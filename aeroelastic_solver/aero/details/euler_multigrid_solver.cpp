#include "euler_multigrid_solver.h"

#include <iostream>
#include <stdexcept>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <vector>

namespace aero {

    void EulerMultigridSolver::coarsen_grid(
        const std::vector<std::vector<double>>& x_fine,
        const std::vector<std::vector<double>>& y_fine,
        std::vector<std::vector<double>>& x_coarse,
        std::vector<std::vector<double>>& y_coarse)
    {
        const int nj_fine = static_cast<int>(x_fine.size());
        const int ni_fine = static_cast<int>(x_fine[0].size());

        // Nodes = cells + 1; both directions need an even cell count so that
        // the last node (index n_cells) is retained and the O-mesh wrap
        // (node 0 == node n_cells) is preserved on the coarse grid.
        if (((ni_fine - 1) % 2) != 0 || ((nj_fine - 1) % 2) != 0) {
            throw std::runtime_error(
                "coarsen_grid: fine cell counts must be even to coarsen.");
        }

        const int ni_coarse = (ni_fine - 1) / 2 + 1;
        const int nj_coarse = (nj_fine - 1) / 2 + 1;

        x_coarse.assign(nj_coarse, std::vector<double>(ni_coarse, 0.0));
        y_coarse.assign(nj_coarse, std::vector<double>(ni_coarse, 0.0));

        for (int J = 0; J < nj_coarse; ++J) {
            for (int I = 0; I < ni_coarse; ++I) {
                x_coarse[J][I] = x_fine[2 * J][2 * I];
                y_coarse[J][I] = y_fine[2 * J][2 * I];
            }
        }
    }

    EulerMultigridSolver::EulerMultigridSolver(
        const std::vector<std::vector<double>>& x_grid,
        const std::vector<std::vector<double>>& y_grid,
        int n_levels)
    {
        if (n_levels < 1) {
            throw std::runtime_error(
                "EulerMultigridSolver: n_levels must be >= 1.");
        }

        levels_.push_back(std::make_unique<EulerSolver>(x_grid, y_grid));

        std::vector<std::vector<double>> cx = x_grid;
        std::vector<std::vector<double>> cy = y_grid;

        for (int k = 1; k < n_levels; ++k) {
            std::vector<std::vector<double>> nx, ny;
            coarsen_grid(cx, cy, nx, ny);
            cx.swap(nx);
            cy.swap(ny);
            levels_.push_back(std::make_unique<EulerSolver>(cx, cy));
        }

        std::cout << "EulerMultigridSolver: built " << levels_.size()
                  << " level(s)\n";
        for (int k = 0; k < num_levels(); ++k) {
            std::cout << "  level " << k
                      << ": " << levels_[k]->ni_cells()
                      << " x " << levels_[k]->nj_cells() << " cells\n";
        }
    }

    void EulerMultigridSolver::set_uniform_state(const Primitive& W) {
        for (auto& lvl : levels_) {
            lvl->set_uniform_state(W);
        }
    }

    void EulerMultigridSolver::restrict_state(
        const EulerSolver& fine, EulerSolver& coarse)
    {
        const int ni_c = coarse.ni_cells();
        const int nj_c = coarse.nj_cells();

        if (fine.ni_cells() != 2 * ni_c || fine.nj_cells() != 2 * nj_c) {
            throw std::runtime_error(
                "restrict_state: coarse grid is not an exact 2x coarsening "
                "of fine grid.");
        }

        for (int J = 0; J < nj_c; ++J) {
            for (int I = 0; I < ni_c; ++I) {

                Conserved sum_U;
                double sum_area = 0.0;

                for (int dj = 0; dj < 2; ++dj) {
                    for (int di = 0; di < 2; ++di) {

                        const int iF = fine.i_start() + 2 * I + di;
                        const int jF = fine.j_start() + 2 * J + dj;

                        const double area = fine.cell_geom(iF, jF).area;

                        sum_area += area;
                        sum_U = sum_U + area * fine.state(iF, jF);
                    }
                }

                coarse.state(coarse.i_start() + I, coarse.j_start() + J) =
                    (1.0 / sum_area) * sum_U;
            }
        }

        coarse.update_periodic_ghost_cells();
    }

    void EulerMultigridSolver::restrict_residual(
        const EulerSolver& fine, EulerSolver& coarse)
    {
        const int ni_c = coarse.ni_cells();
        const int nj_c = coarse.nj_cells();

        if (fine.ni_cells() != 2 * ni_c || fine.nj_cells() != 2 * nj_c) {
            throw std::runtime_error(
                "restrict_residual: coarse grid is not an exact 2x coarsening "
                "of fine grid.");
        }

        for (int J = 0; J < nj_c; ++J) {
            for (int I = 0; I < ni_c; ++I) {

                // Sum the integrated (area * per-volume) fine residuals over
                // the 2x2 block. Internal face fluxes cancel, so this is the
                // integrated flux balance of the coarse cell.
                Conserved sum_R;

                for (int dj = 0; dj < 2; ++dj) {
                    for (int di = 0; di < 2; ++di) {

                        const int iF = fine.i_start() + 2 * I + di;
                        const int jF = fine.j_start() + 2 * J + dj;

                        const double area = fine.cell_geom(iF, jF).area;
                        sum_R = sum_R + area * fine.residual(iF, jF);
                    }
                }

                // Store back per (coarse) unit volume, dividing by the coarse
                // cell's own area so area_coarse * R_coarse == sum_R exactly.
                const int iC = coarse.i_start() + I;
                const int jC = coarse.j_start() + J;
                const double area_coarse = coarse.cell_geom(iC, jC).area;

                coarse.residual(iC, jC) = (1.0 / area_coarse) * sum_R;
            }
        }
    }

    void EulerMultigridSolver::compute_forcing(
        EulerSolver& fine, EulerSolver& coarse,
        const Primitive& Winf,
        const JSTParameters& fine_parameters,
        const JSTParameters& coarse_parameters)
    {
        // 1. Fine-level FAS defect d_h = R_h(u_h) - forcing_h, written back into
        //    fine.R_ so restrict_residual transfers the defect. On the finest
        //    level forcing_h == 0 (d_h == R_h); the subtraction only matters on
        //    intermediate levels once recursion is added (step 12). The fine
        //    residual uses the fine operator (full JST on the finest level).
        fine.compute_full_residual_jst(Winf, fine_parameters);
        {
            std::vector<Conserved>&       Rf = fine.residual_field();
            const std::vector<Conserved>& Ff = fine.forcing_field();
            for (int j = fine.j_start(); j < fine.j_end(); ++j) {
                for (int i = fine.i_start(); i < fine.i_end(); ++i) {
                    const int id = fine.cell_index(i, j);
                    Rf[id] = Rf[id] - Ff[id];
                }
            }
        }

        // 2. Restrict the fine state: coarse.U = I u_h. (Updates coarse ghosts.)
        restrict_state(fine, coarse);

        // 3. Coarse operator on the restricted state: coarse.R_ = R_H(I u_h),
        //    using the coarse operator (first-order on coarse levels). This is
        //    the operator the coarse smoother will drive R_H - tau to zero
        //    against, so it must match smooth_level's operator. Stash R_H(I u_h)
        //    in coarse.forcing_ before restrict_residual overwrites R_.
        coarse.compute_full_residual_jst(Winf, coarse_parameters);
        {
            const std::vector<Conserved>& Rc = coarse.residual_field();
            std::vector<Conserved>&       Fc = coarse.forcing_field();
            for (int j = coarse.j_start(); j < coarse.j_end(); ++j) {
                for (int i = coarse.i_start(); i < coarse.i_end(); ++i) {
                    const int id = coarse.cell_index(i, j);
                    Fc[id] = Rc[id];
                }
            }
        }

        // 4. Restrict the fine defect: coarse.R_ = I d_h.
        restrict_residual(fine, coarse);

        // 5. tau = R_H(I u_h) - I d_h, left in coarse.forcing_.
        {
            const std::vector<Conserved>& Rc = coarse.residual_field();
            std::vector<Conserved>&       Fc = coarse.forcing_field();
            for (int j = coarse.j_start(); j < coarse.j_end(); ++j) {
                for (int i = coarse.i_start(); i < coarse.i_end(); ++i) {
                    const int id = coarse.cell_index(i, j);
                    Fc[id] = Fc[id] - Rc[id];
                }
            }
        }
    }

    void EulerMultigridSolver::smooth_level(
        int lvl, const Primitive& Winf, const JSTParameters& parameters,
        double cfl, int n_sweeps)
    {
        EulerSolver& L = level(lvl);
        for (int s = 0; s < n_sweeps; ++s) {
            // Freeze D_ for this sweep's five RK stages (recompute_dissipation
            // defaults to true), then take one RK5 + IRS step. The FAS forcing
            // (tau) is subtracted inside each stage; it is zero on the finest
            // level. Matches the single-grid driver's per-iteration sequence.
            L.compute_full_residual_jst(Winf, parameters);
            L.rk5_pseudo_step(Winf, parameters, cfl);
        }
    }

    void EulerMultigridSolver::prolong_correction(
        EulerSolver& fine, const EulerSolver& coarse,
        const std::vector<Conserved>& coarse_pre, double omega, bool bilinear)
    {
        const int ni_c = coarse.ni_cells();
        const int nj_c = coarse.nj_cells();

        if (fine.ni_cells() != 2 * ni_c || fine.nj_cells() != 2 * nj_c) {
            throw std::runtime_error(
                "prolong_correction: coarse grid is not an exact 2x coarsening "
                "of fine grid.");
        }

        const std::vector<Conserved>& Uc = coarse.state_field();
        std::vector<Conserved>&       Uf = fine.state_field();

        // FAS correction c_H(I,J) = (smoothed coarse state) - (restricted
        // pre-smoothing state I u_h), at a physical coarse cell. i is periodic
        // (wrap the index), j is wall/farfield (clamp = zero-gradient).
        auto cH_at = [&](int I, int J) -> Conserved {
            const int Ip = ((I % ni_c) + ni_c) % ni_c;
            const int Jp = std::max(0, std::min(nj_c - 1, J));
            const int id = coarse.cell_index(coarse.i_start() + Ip,
                                             coarse.j_start() + Jp);
            return Uc[id] - coarse_pre[id];
        };

        for (int J = 0; J < nj_c; ++J) {
            for (int I = 0; I < ni_c; ++I) {
                for (int dj = 0; dj < 2; ++dj) {
                    for (int di = 0; di < 2; ++di) {

                        Conserved c;
                        if (bilinear) {
                            // Cell-centered bilinear: the fine cell is offset a
                            // quarter coarse-cell toward (I_n, J_n); interpolate
                            // from the 2x2 bracketing coarse cells with the
                            // ideal-grid weights 9/16, 3/16, 3/16, 1/16.
                            const int I_n = (di == 0) ? I - 1 : I + 1;
                            const int J_n = (dj == 0) ? J - 1 : J + 1;
                            c =   (9.0 / 16.0) * cH_at(I,   J)
                                + (3.0 / 16.0) * cH_at(I_n, J)
                                + (3.0 / 16.0) * cH_at(I,   J_n)
                                + (1.0 / 16.0) * cH_at(I_n, J_n);
                        }
                        else {
                            // Piecewise injection: all four fine cells take the
                            // same coarse correction.
                            c = cH_at(I, J);
                        }

                        const int iF  = fine.i_start() + 2 * I + di;
                        const int jF  = fine.j_start() + 2 * J + dj;
                        const int idF = fine.cell_index(iF, jF);
                        Uf[idF] = Uf[idF] + omega * c;
                    }
                }
            }
        }

        fine.update_periodic_ghost_cells();
    }

    void EulerMultigridSolver::prolong_solution(
        EulerSolver& fine, const EulerSolver& coarse, bool bilinear)
    {
        const int ni_c = coarse.ni_cells();
        const int nj_c = coarse.nj_cells();

        if (fine.ni_cells() != 2 * ni_c || fine.nj_cells() != 2 * nj_c) {
            throw std::runtime_error(
                "prolong_solution: coarse grid is not an exact 2x coarsening "
                "of fine grid.");
        }

        const std::vector<Conserved>& Uc = coarse.state_field();
        std::vector<Conserved>&       Uf = fine.state_field();

        // Coarse solution at physical cell (I,J); i periodic (wrap), j clamped
        // (zero-gradient at wall/farfield). Same access pattern as the bilinear
        // correction stencil, but on the absolute state.
        auto uH_at = [&](int I, int J) -> Conserved {
            const int Ip = ((I % ni_c) + ni_c) % ni_c;
            const int Jp = std::max(0, std::min(nj_c - 1, J));
            const int id = coarse.cell_index(coarse.i_start() + Ip,
                                             coarse.j_start() + Jp);
            return Uc[id];
        };

        for (int J = 0; J < nj_c; ++J) {
            for (int I = 0; I < ni_c; ++I) {
                for (int dj = 0; dj < 2; ++dj) {
                    for (int di = 0; di < 2; ++di) {

                        Conserved u;
                        if (bilinear) {
                            const int I_n = (di == 0) ? I - 1 : I + 1;
                            const int J_n = (dj == 0) ? J - 1 : J + 1;
                            u =   (9.0 / 16.0) * uH_at(I,   J)
                                + (3.0 / 16.0) * uH_at(I_n, J)
                                + (3.0 / 16.0) * uH_at(I,   J_n)
                                + (1.0 / 16.0) * uH_at(I_n, J_n);
                        }
                        else {
                            u = uH_at(I, J);
                        }

                        const int iF  = fine.i_start() + 2 * I + di;
                        const int jF  = fine.j_start() + 2 * J + dj;
                        Uf[fine.cell_index(iF, jF)] = u;
                    }
                }
            }
        }

        fine.update_periodic_ghost_cells();
    }

    void EulerMultigridSolver::v_cycle(
        int lvl, int finest_lvl, const Primitive& Winf,
        const JSTParameters& parameters, const JSTParameters& coarse_parameters,
        double cfl_fine, double coarse_cfl_decay,
        int n_pre, int n_coarse, int n_post,
        double correction_omega, double omega_decay, bool bilinear,
        int cycle_gamma)
    {
        // Operator for this level: full JST on this cycle's finest level (the
        // real grid being solved -- level 0 in a normal run, or the temporary
        // FMG-ramp top level), first-order on every coarser level below it.
        const JSTParameters& p_lvl =
            (lvl == finest_lvl) ? parameters : coarse_parameters;

        // Per-level throttle: full CFL / omega at this cycle's finest level,
        // geometrically reduced for each coarsening below it. Measured relative
        // to finest_lvl so an FMG ramp top level always starts at full CFL.
        const int    depth     = lvl - finest_lvl;
        const double cfl_lvl   = cfl_fine         * std::pow(coarse_cfl_decay, depth);
        const double omega_lvl = correction_omega * std::pow(omega_decay,      depth);

        // Base case: coarsest level, approximately solved by smoothing.
        if (lvl == num_levels() - 1) {
            smooth_level(lvl, Winf, p_lvl, cfl_lvl, n_coarse);
            return;
        }

        // 1. Pre-smooth this level.
        smooth_level(lvl, Winf, p_lvl, cfl_lvl, n_pre);

        // 2. Restrict state and build tau on lvl+1. fine_params = this level's
        //    operator (so the restricted defect R_lvl - forcing_lvl is formed
        //    with it); coarse_params = first-order (lvl+1 is always a coarse
        //    level).
        compute_forcing(level(lvl), level(lvl + 1), Winf, p_lvl, coarse_parameters);

        // 3. Save the restricted pre-smoothing state I u_lvl.
        std::vector<Conserved> coarse_pre = level(lvl + 1).state_field();

        // 4. Recurse on the coarser level cycle_gamma times (gamma = 1 is the
        //    V-cycle; gamma = 2 is the W-cycle). The coarse forcing (tau) was
        //    fixed in step 2 and is held across the repeats -- each repeat just
        //    relaxes the coarse problem further toward R_{lvl+1} = tau. Because
        //    coarse_pre was captured before the first visit (step 3), the
        //    correction prolonged in step 5 is the NET change over all gamma
        //    visits. gamma only deepens the recursion below this level, so it
        //    has no effect with only 2 levels (lvl+1 is the coarsest base case).
        for (int g = 0; g < cycle_gamma; ++g) {
            v_cycle(lvl + 1, finest_lvl, Winf, parameters, coarse_parameters,
                    cfl_fine, coarse_cfl_decay, n_pre, n_coarse, n_post,
                    correction_omega, omega_decay, bilinear, cycle_gamma);
        }

        // 5. Prolong the (per-level under-relaxed) correction back to this level.
        prolong_correction(level(lvl), level(lvl + 1), coarse_pre,
                           omega_lvl, bilinear);

        // 6. Post-smooth this level.
        smooth_level(lvl, Winf, p_lvl, cfl_lvl, n_post);
    }

    int EulerMultigridSolver::run_v_cycles(
        const Primitive& Winf, const JSTParameters& parameters,
        double cfl, int max_cycles, int print_every, double target_residual,
        int n_pre, int n_coarse, int n_post, double correction_omega,
        double omega_decay, double coarse_cfl_decay, bool bilinear_prolong,
        int cycle_gamma)
    {
        if (num_levels() < 2) {
            throw std::runtime_error(
                "run_v_cycles requires at least 2 grid levels.");
        }

        // The coarse level runs a robust first-order dissipation operator
        // (constant epsilon_2, no 4th-difference, no pressure-sensor switching)
        // to stop the nonlinear JST coarse correction from overshooting into a
        // non-physical state during the transonic cold-start transient. The
        // fine level keeps the full (second-order) JST scheme in `parameters`.
        JSTParameters coarse_parameters = parameters;
        coarse_parameters.first_order = true;

        const char* cycle_name =
            (cycle_gamma <= 1) ? "V" : (cycle_gamma == 2 ? "W" : "mu");
        std::cout << "\nRunning " << num_levels()
                  << "-level FAS " << cycle_name << "-cycles (RK5 + JST + IRS)"
                  << "  [gamma=" << cycle_gamma << "]\n";
        std::cout << "---------------------------------------------\n";
        std::cout << "  max_cycles      = " << max_cycles << "\n";
        std::cout << "  target_residual = " << std::scientific << target_residual << "\n";
        std::cout << "  CFL             = " << cfl << "\n";
        std::cout << "  smoothing       = " << n_pre << " pre / "
                  << n_coarse << " coarse / " << n_post << " post\n";
        std::cout << "  coarse dissip.  = first-order (eps2=0.5, eps4=0)\n";
        std::cout << "  correction omega= " << correction_omega
                  << " (decay " << omega_decay << "/level)\n";
        std::cout << "  coarse cfl decay= " << coarse_cfl_decay << "/level\n";
        std::cout << "  prolongation    = "
                  << (bilinear_prolong ? "bilinear" : "piecewise") << "\n";
        std::cout << "  per-level schedule (cfl, omega):\n";
        for (int k = 0; k < num_levels(); ++k) {
            std::cout << "    level " << k
                      << ": cfl=" << (cfl * std::pow(coarse_cfl_decay, k))
                      << "  omega=" << (correction_omega * std::pow(omega_decay, k))
                      << "\n";
        }
        std::cout << "\n";

        const double alpha_rad = std::atan2(Winf.v, Winf.u);
        const auto start_time = std::chrono::steady_clock::now();

        int  final_cycle = 0;
        bool converged   = false;
        long fine_sweeps = 0;

        for (int cyc = 0; cyc < max_cycles; ++cyc) {

            // True (unsmoothed) fine residual for monitoring and convergence.
            fine().compute_full_residual_jst(Winf, parameters);
            const double res = fine().residual_linf_current();

            if (cyc == 0 || cyc % print_every == 0 || cyc == max_cycles - 1) {
                const auto now = std::chrono::steady_clock::now();
                const double mins =
                    std::chrono::duration<double>(now - start_time).count() / 60.0;
                const AeroCoefficients c =
                    fine().compute_aero_coefficients(Winf, alpha_rad);

                std::cout << std::scientific << std::setprecision(6);
                std::cout << "  cycle " << std::setw(5) << cyc
                          << "  residual_linf = " << res
                          << "  CL = " << c.cl
                          << "  CD = " << c.cd
                          << "  CM = " << c.cm
                          << "  fine_sweeps = " << fine_sweeps
                          << "  elapsed = " << std::fixed << std::setprecision(2)
                          << mins << " min"
                          << std::scientific << std::setprecision(6) << "\n";
            }

            if (res <= target_residual) {
                converged   = true;
                final_cycle = cyc;
                break;
            }

            // CFL ramp over the first n_ramp cycles to absorb the cold-start
            // transient, mirroring the single-grid driver's per-iteration ramp.
            const int    n_ramp    = 50;
            const double cfl_start = 0.3;
            const double frac =
                (cyc >= n_ramp) ? 1.0 : static_cast<double>(cyc) / n_ramp;
            const double cfl_now =
                std::min(cfl, cfl_start + frac * (cfl - cfl_start));

            // One recursive FAS V-cycle over the whole hierarchy. The per-level
            // CFL/omega decay, first-order coarse dissipation, correction
            // under-relaxation, and prolongation operator are applied at every
            // level inside v_cycle().
            v_cycle(0, /*finest_lvl=*/0, Winf, parameters, coarse_parameters,
                    cfl_now, coarse_cfl_decay, n_pre, n_coarse, n_post,
                    correction_omega, omega_decay, bilinear_prolong, cycle_gamma);

            // Count only level-0 (fine) smoothing sweeps as the work metric;
            // coarser levels are geometrically cheaper.
            fine_sweeps += n_pre + n_post;
            final_cycle  = cyc;
        }

        fine().compute_full_residual_jst(Winf, parameters);
        const double final_res = fine().residual_linf_current();
        const auto   end_time  = std::chrono::steady_clock::now();
        const double total_mins =
            std::chrono::duration<double>(end_time - start_time).count() / 60.0;

        std::cout << "  final residual_linf = "
                  << std::scientific << std::setprecision(6) << final_res << "\n";
        std::cout << "  final cycle         = " << final_cycle << "\n";
        std::cout << "  fine sweeps         = " << fine_sweeps << "\n";
        std::cout << "  total runtime       = "
                  << std::fixed << std::setprecision(2) << total_mins << " min\n";
        std::cout << "  status              = "
                  << (converged ? "converged" : "NOT converged") << "\n";

        return converged ? final_cycle : max_cycles;
    }

    void EulerMultigridSolver::fmg_ramp(
        const Primitive& Winf, const JSTParameters& parameters,
        const JSTParameters& coarse_parameters, double cfl, int n_fmg_cycles,
        int n_pre, int n_coarse, int n_post, double correction_omega,
        double omega_decay, double coarse_cfl_decay, bool bilinear,
        int cycle_gamma)
    {
        const auto ramp_start = std::chrono::steady_clock::now();

        // Work from the coarsest level up to (but not including) the finest:
        // converge level L (as the temporary finest grid, so it uses the full
        // operator -- a first-order ramp solution makes too poor an initial
        // guess), then prolong its solution to level L-1. After the loop,
        // level 0 holds the FMG initial guess.
        for (int L = num_levels() - 1; L >= 1; --L) {

            // Each stage starts from a freshly-prolonged (coarser-grid) field,
            // so the finer grid sharpens features -- especially a transonic
            // shock -- over the first cycles. Jumping straight to full CFL on a
            // fine grid overshoots and can drive a cell non-physical (the cause
            // of the 512^2 transonic crash). So ramp the CFL within EVERY stage,
            // exactly as run_v_cycles ramps the main solve.
            const double cfl_start = std::min(0.3, cfl);
            const int    n_ramp    = 20;
            auto stage_cfl = [&](int s, double cfl_target) {
                const double frac =
                    (s >= n_ramp) ? 1.0 : static_cast<double>(s) / n_ramp;
                return std::min(cfl_target,
                                cfl_start + frac * (cfl_target - cfl_start));
            };

            if (L == num_levels() - 1) {
                // Coarsest level: no coarser grid to recurse to, so just smooth.
                // This is FMG's only cold start (from uniform freestream).
                // Deliberately use the robust FIRST-ORDER operator here, not the
                // full JST one: on the coarsest grid the 2nd-order accuracy is
                // negligible (its solution is immediately refined by the
                // full-JST stage one level up), while full JST's pressure sensor
                // can throw on a non-physical coarse-cell pressure. The constant
                // eps2=0.5 first-order operator has a TIGHTER CFL limit than full
                // JST, so target the throttled CFL this level uses in the main
                // V-cycle (cfl * decay^depth), not the full fine-grid CFL.
                const double cfl_target =
                    cfl * std::pow(coarse_cfl_decay, num_levels() - 1);
                for (int s = 0; s < n_fmg_cycles; ++s) {
                    smooth_level(L, Winf, coarse_parameters,
                                 stage_cfl(s, cfl_target), 1);
                }
            }
            else {
                // Intermediate level: run FAS V-cycles with L as the temporary
                // finest grid (full operator at L, first-order / decayed below
                // it), ramping the CFL toward the full value over the first
                // cycles of the stage.
                for (int c = 0; c < n_fmg_cycles; ++c) {
                    v_cycle(L, /*finest_lvl=*/L, Winf, parameters,
                            coarse_parameters, stage_cfl(c, cfl), coarse_cfl_decay,
                            n_pre, n_coarse, n_post, correction_omega,
                            omega_decay, bilinear, cycle_gamma);
                }
            }

            // Report this stage's result (and reveal the crash level, if any,
            // since the ramp is otherwise silent). Use the operator the stage
            // used so the residual is meaningful.
            const JSTParameters& p_used =
                (L == num_levels() - 1) ? coarse_parameters : parameters;
            level(L).compute_full_residual_jst(Winf, p_used);
            const double stage_res = level(L).residual_linf_current();
            const double mins =
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - ramp_start).count() / 60.0;
            std::cout << std::scientific << std::setprecision(6)
                      << "  FMG level " << L << " ("
                      << level(L).ni_cells() << "x" << level(L).nj_cells()
                      << "): residual_linf = " << stage_res
                      << std::fixed << std::setprecision(2)
                      << "  (" << mins << " min)\n" << std::flush;

            // Seed the next-finer level with this level's converged solution.
            prolong_solution(level(L - 1), level(L), bilinear);
        }
    }

    int EulerMultigridSolver::run_full_multigrid(
        const Primitive& Winf, const JSTParameters& parameters,
        double cfl, int max_cycles, int print_every, double target_residual,
        int n_pre, int n_coarse, int n_post, int n_fmg_cycles,
        double correction_omega, double omega_decay, double coarse_cfl_decay,
        bool bilinear_prolong, int cycle_gamma)
    {
        if (num_levels() < 2) {
            throw std::runtime_error(
                "run_full_multigrid requires at least 2 grid levels.");
        }

        JSTParameters coarse_parameters = parameters;
        coarse_parameters.first_order = true;

        std::cout << "\nFull Multigrid (FMG) start-up\n";
        std::cout << "---------------------------------------------\n";
        std::cout << "  ramp levels     = coarsest (level "
                  << (num_levels() - 1) << ") -> finest (level 0)\n";
        std::cout << "  cycles per level= " << n_fmg_cycles << "\n";
        std::cout << "  prolongation    = "
                  << (bilinear_prolong ? "bilinear" : "piecewise")
                  << " (full solution transfer)\n\n";

        const auto fmg_start = std::chrono::steady_clock::now();

        fmg_ramp(Winf, parameters, coarse_parameters, cfl, n_fmg_cycles,
                 n_pre, n_coarse, n_post, correction_omega, omega_decay,
                 coarse_cfl_decay, bilinear_prolong, cycle_gamma);

        fine().compute_full_residual_jst(Winf, parameters);
        const double res_after_fmg = fine().residual_linf_current();
        const auto   fmg_end = std::chrono::steady_clock::now();
        const double fmg_mins =
            std::chrono::duration<double>(fmg_end - fmg_start).count() / 60.0;

        std::cout << std::scientific << std::setprecision(6)
                  << "  FMG ramp done: fine residual_linf = " << res_after_fmg
                  << std::fixed << std::setprecision(2)
                  << "  (" << fmg_mins << " min)\n";

        // Level 0 now holds the FMG initial guess; drive it to convergence.
        return run_v_cycles(Winf, parameters, cfl, max_cycles, print_every,
                            target_residual, n_pre, n_coarse, n_post,
                            correction_omega, omega_decay, coarse_cfl_decay,
                            bilinear_prolong, cycle_gamma);
    }
}
