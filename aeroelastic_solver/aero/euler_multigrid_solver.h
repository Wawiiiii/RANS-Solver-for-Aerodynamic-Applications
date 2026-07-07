#pragma once

// FAS multigrid engine for the Euler kernel. Owns a hierarchy of EulerSolver
// levels (0 = finest) built by every-other-node coarsening of the O-mesh, so
// coarse nodes are an exact subset of fine nodes (properly nested grids -> clean
// transfer operators). Drives V/W-cycles and Full Multigrid start-up.

#include "euler_solver.h"

#include <memory>
#include <vector>

namespace aero {

    class EulerMultigridSolver {
    public:
        // Build n_levels by repeated 2x coarsening of the given finest node grid.
        EulerMultigridSolver(
            const std::vector<std::vector<double>>& x_grid,
            const std::vector<std::vector<double>>& y_grid,
            int n_levels
        );

        int num_levels() const { return static_cast<int>(levels_.size()); }

        EulerSolver&       fine()       { return *levels_.front(); }
        const EulerSolver& fine() const { return *levels_.front(); }

        // Set the same uniform state on every level.
        void set_uniform_state(const Primitive& W);

        // Drive the finest level with FAS cycles until residual_linf reaches
        // target_residual or max_cycles. Each cycle: n_pre fine smooths -> build
        // tau + recurse to coarser levels -> prolong correction -> n_post smooths.
        // CFL is ramped over the first cycles to absorb the cold start.
        // Per-level throttle (high near fine, decayed toward coarse, where larger
        // cells make corrections riskier):
        //   cfl_lvl   = cfl              * coarse_cfl_decay^lvl
        //   omega_lvl = correction_omega * omega_decay^lvl
        // cycle_gamma = FAS recursion factor: 1 = V-cycle, 2 = W-cycle (visit
        // each coarser level twice; the transonic speed lever). Only affects 3+
        // levels. Returns the convergence cycle (or max_cycles if not converged).
        int run_v_cycles(
            const Primitive& Winf,
            const JSTParameters& parameters,
            double cfl,
            int max_cycles,
            int print_every,
            double target_residual,
            int n_pre,
            int n_coarse,
            int n_post,
            double correction_omega = 1.0,
            double omega_decay = 1.0,
            double coarse_cfl_decay = 1.0,
            bool bilinear_prolong = true,
            int cycle_gamma = 1
        );

        // Full Multigrid: converge coarsest -> prolong SOLUTION up as the next
        // finer level's initial guess -> repeat to finest, then run_v_cycles to
        // target. Removes the bulk of the cold-start transient. n_fmg_cycles =
        // work per ramp level. Only changes the initial guess, not the answer.
        // Remaining parameters carry the run_v_cycles meaning. Returns the
        // run_v_cycles convergence cycle.
        int run_full_multigrid(
            const Primitive& Winf,
            const JSTParameters& parameters,
            double cfl,
            int max_cycles,
            int print_every,
            double target_residual,
            int n_pre,
            int n_coarse,
            int n_post,
            int n_fmg_cycles,
            double correction_omega = 1.0,
            double omega_decay = 1.0,
            double coarse_cfl_decay = 1.0,
            bool bilinear_prolong = true,
            int cycle_gamma = 1
        );

    private:
        EulerSolver&       level(int k)       { return *levels_.at(k); }
        const EulerSolver& level(int k) const { return *levels_.at(k); }

        // Coarsen a node grid by taking every other node in both directions.
        static void coarsen_grid(
            const std::vector<std::vector<double>>& x_fine,
            const std::vector<std::vector<double>>& y_fine,
            std::vector<std::vector<double>>& x_coarse,
            std::vector<std::vector<double>>& y_coarse
        );

        // Restrict state: fine-area-weighted average over each 2x2 fine block
        // (preserves the area-weighted mean of U). Updates coarse ghosts.
        static void restrict_state(const EulerSolver& fine, EulerSolver& coarse);

        // Restrict residual conservatively: area_coarse*R_coarse = sum of the
        // 2x2 fine integrated residuals. This is the I*R_h term of tau.
        static void restrict_residual(const EulerSolver& fine, EulerSolver& coarse);

        // Build FAS forcing tau = R_H(I u_h) - I(R_h - forcing_h) into
        // coarse.forcing_field(). fine_parameters drives R_h; coarse_parameters
        // (first_order) drives R_H and must match the coarse smoother's operator.
        static void compute_forcing(
            EulerSolver& fine,
            EulerSolver& coarse,
            const Primitive& Winf,
            const JSTParameters& fine_parameters,
            const JSTParameters& coarse_parameters
        );

        // n_sweeps RK5 + IRS relaxations at fixed cfl (no startup ramp). Each
        // sweep freezes D_ then takes one rk5_pseudo_step (subtracts tau).
        void smooth_level(
            int lvl,
            const Primitive& Winf,
            const JSTParameters& parameters,
            double cfl,
            int n_sweeps
        );

        // Prolong the FAS correction c_H = coarse.U - coarse_pre and add it:
        // u_h += omega * P c_H. omega < 1 under-relaxes (bounds a large coarse
        // excursion -> transonic stability). bilinear: true = cell-centered
        // bilinear (9/16,3/16,3/16,1/16, periodic i, clamped j); false = injection.
        static void prolong_correction(
            EulerSolver& fine,
            const EulerSolver& coarse,
            const std::vector<Conserved>& coarse_pre,
            double omega = 1.0,
            bool bilinear = true
        );

        // Prolong the full coarse SOLUTION, overwriting the fine state
        // (u_h <- P u_H). Same operator P as prolong_correction. Used by FMG.
        static void prolong_solution(
            EulerSolver& fine,
            const EulerSolver& coarse,
            bool bilinear = true
        );

        // FMG ramp: coarsest -> finest, prolonging the solution up each step.
        // Leaves level 0 holding the FMG initial guess. Coarsest stage smooths
        // (first-order, internal CFL ramp); intermediate stages run FAS cycles.
        void fmg_ramp(
            const Primitive& Winf,
            const JSTParameters& parameters,
            const JSTParameters& coarse_parameters,
            double cfl,
            int n_fmg_cycles,
            int n_pre,
            int n_coarse,
            int n_post,
            double correction_omega,
            double omega_decay,
            double coarse_cfl_decay,
            bool bilinear,
            int cycle_gamma
        );

        // One recursive FAS cycle from level lvl. finest_lvl = this cycle's real
        // top grid (full JST there, first_order below; throttle measured relative
        // to it). cycle_gamma = how many times the coarser level is visited (tau
        // held fixed across visits; prolonged correction = net change).
        void v_cycle(
            int lvl,
            int finest_lvl,
            const Primitive& Winf,
            const JSTParameters& parameters,
            const JSTParameters& coarse_parameters,
            double cfl_fine,
            double coarse_cfl_decay,
            int n_pre,
            int n_coarse,
            int n_post,
            double correction_omega,
            double omega_decay,
            bool bilinear,
            int cycle_gamma
        );

        std::vector<std::unique_ptr<EulerSolver>> levels_;
    };
}
