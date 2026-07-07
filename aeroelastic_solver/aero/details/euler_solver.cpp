#include "euler_solver.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace aero {

    namespace {

        // Non-periodic tridiagonal solve (Thomas algorithm). Coefficients
        // sub/diag/sup have length n (sub[0] and sup[n-1] are unused). rhs is
        // overwritten with the solution; cprime is scratch of length >= n.
        void thomas_scalar(
            const std::vector<double>& sub,
            const std::vector<double>& diag,
            const std::vector<double>& sup,
            std::vector<double>& rhs,
            int n,
            std::vector<double>& cprime)
        {
            cprime[0] = sup[0] / diag[0];
            rhs[0]    = rhs[0] / diag[0];

            for (int k = 1; k < n; ++k) {
                const double m = 1.0 / (diag[k] - sub[k] * cprime[k - 1]);
                cprime[k] = sup[k] * m;
                rhs[k]    = (rhs[k] - sub[k] * rhs[k - 1]) * m;
            }

            for (int k = n - 2; k >= 0; --k) {
                rhs[k] -= cprime[k] * rhs[k + 1];
            }
        }

        // Periodic (cyclic) tridiagonal solve via Sherman-Morrison.
        // Row k: sub[k] x[k-1] + diag[k] x[k] + sup[k] x[k+1] = rhs[k] (mod n).
        // beta  = sub[0]    couples x[n-1] into row 0   (top-right corner).
        // alpha = sup[n-1]  couples x[0]   into row n-1 (bottom-left corner).
        // rhs is overwritten with the solution; bb/u/cprime are scratch.
        void cyclic_scalar(
            const std::vector<double>& sub,
            const std::vector<double>& diag,
            const std::vector<double>& sup,
            std::vector<double>& rhs,
            int n,
            std::vector<double>& bb,
            std::vector<double>& u,
            std::vector<double>& cprime)
        {
            if (n == 1) {
                rhs[0] = rhs[0] / diag[0];
                return;
            }

            const double beta  = sub[0];
            const double alpha = sup[n - 1];
            const double gamma = -diag[0];

            for (int k = 0; k < n; ++k) bb[k] = diag[k];
            bb[0]     = diag[0] - gamma;
            bb[n - 1] = diag[n - 1] - alpha * beta / gamma;

            // Solve A' x = rhs (A' = tridiag with modified diagonal bb).
            thomas_scalar(sub, bb, sup, rhs, n, cprime);

            // Solve A' z = u, with u = gamma e_0 + alpha e_{n-1}.
            for (int k = 0; k < n; ++k) u[k] = 0.0;
            u[0]     = gamma;
            u[n - 1] = alpha;
            thomas_scalar(sub, bb, sup, u, n, cprime);

            const double fact =
                (rhs[0] + beta * rhs[n - 1] / gamma) /
                (1.0 + u[0] + beta * u[n - 1] / gamma);

            for (int k = 0; k < n; ++k) rhs[k] -= fact * u[k];
        }
    }

    // Internal flux/boundary helpers (file-local; not exposed in the header).
    double    total_energy_per_mass(const Primitive& W);
    Conserved normal_flux(const Conserved& U, double nx, double ny);
    Conserved central_flux(const Conserved& UL, const Conserved& UR, double nx, double ny);
    double    face_spectral_radius(const Conserved& UL, const Conserved& UR, double nx, double ny);
    double    pressure_sensor(double p_minus, double p_center, double p_plus);
    Conserved averaged_state_flux(const Primitive& WL, const Primitive& WR, double nx, double ny);
    Conserved jst_artificial_flux(
        const Conserved& U_minus, const Conserved& U_left,
        const Conserved& U_right, const Conserved& U_plus,
        double lambda_face, double epsilon_2, double epsilon_4);
    Primitive reflected_wall_state(const Primitive& Wi, double nx, double ny);
    Primitive farfield_riemann_state(const Primitive& Wi, const Primitive& Winf, double nx, double ny);

    bool is_physical(const Primitive& W) {
        return W.rho > 0.0 && W.p > 0.0;
    }

    double total_energy_per_mass(const Primitive& W) { 

        if (!is_physical(W)) { 
            throw std::runtime_error("Non-Physical primitive state in total_energy_per_mass.");
        }

        const double kinetic = 0.5 * (W.u * W.u + W.v * W.v); 
        const double internal = W.p / ((GAMMA - 1.0) * W.rho);

        return internal + kinetic;
    }

    Conserved primitive_to_conserved(const Primitive& W) {

        if (!is_physical(W))
        {
            throw std::runtime_error("Non-physical primitive state in primitive_to_conserved.");
        }

        const double E = total_energy_per_mass(W);

        Conserved U;
        U.rho  = W.rho;
        U.rhou = W.rho * W.u;
        U.rhov = W.rho * W.v;
        U.rhoE = W.rho * E;

        return U;
    }

    Primitive conserved_to_primitive(const Conserved& U) {
        if (U.rho <= 0.0) {
            throw std::runtime_error("Non-physical density in conserved_to_primitive.");
        }

        Primitive W;

        W.rho = U.rho;
        W.u   = U.rhou / U.rho;
        W.v   = U.rhov / U.rho;

        const double E = U.rhoE / U.rho;
        const double kinetic = 0.5 * (W.u * W.u + W.v * W.v);

        W.p = (GAMMA - 1.0) * U.rho * (E - kinetic);

        if (!is_physical(W)) {
            throw std::runtime_error("Non-physical pressure in conserved_to_primitive.");
        }

        return W;
    }

    double sound_speed(const Primitive& W) {

        if (!is_physical(W)) {
            throw std::runtime_error("Non-physical primitive state in sound_speed.");
        }

        return std::sqrt(GAMMA * W.p / W.rho);
    }

    Conserved operator+(const Conserved& a, const Conserved& b) {

        return {
            a.rho  + b.rho,
            a.rhou + b.rhou,
            a.rhov + b.rhov,
            a.rhoE + b.rhoE
        };
    }

    Conserved operator-(const Conserved& a, const Conserved& b) {
        
        return {
            a.rho  - b.rho,
            a.rhou - b.rhou,
            a.rhov - b.rhov,
            a.rhoE - b.rhoE
        };
    }

    Conserved operator*(double s, const Conserved& a) {

        return {
            s * a.rho,
            s * a.rhou,
            s * a.rhov,
            s * a.rhoE
        };
    }

    Conserved operator*(const Conserved& a, double s) {
        return s * a;
    }

    Conserved normal_flux(const Conserved& U, double nx, double ny) {
        const Primitive W = conserved_to_primitive(U);

        const double rho = W.rho;
        const double u   = W.u;
        const double v   = W.v;
        const double p   = W.p;

        const double un = u * nx + v * ny;

        Conserved F;

        F.rho  = rho * un;
        F.rhou = rho * u * un + p * nx;
        F.rhov = rho * v * un + p * ny;
        F.rhoE = (U.rhoE + p) * un;

        return F;
    }
        EulerSolver::EulerSolver(const std::vector<std::vector<double>>& x_grid, const std::vector<std::vector<double>>& y_grid) {
        validate_grid(x_grid, y_grid);

        nj_nodes_ = static_cast<int>(x_grid.size());
        ni_nodes_ = static_cast<int>(x_grid[0].size());

        ni_cells_ = ni_nodes_ - 1;
        nj_cells_ = nj_nodes_ - 1;

        ni_total_ = ni_cells_ + 2 * ng_;
        nj_total_ = nj_cells_ + 2 * ng_;

        i_start_ = ng_;
        i_end_   = ng_ + ni_cells_;

        j_start_ = ng_;
        j_end_   = ng_ + nj_cells_;

        geom_.resize(ni_total_ * nj_total_);
        U_.resize(ni_total_ * nj_total_); 
        R_.resize(ni_total_ * nj_total_);
        D_.resize(ni_total_ * nj_total_);
        forcing_.resize(ni_total_ * nj_total_);
        eps_i_.resize(ni_total_ * nj_total_, 0.0);
        eps_j_.resize(ni_total_ * nj_total_, 0.0);

        build_geometry(x_grid, y_grid);

        check_geometry();
    }

    Conserved averaged_state_flux( const Primitive& WL, const Primitive& WR, double nx, double ny) {
        Primitive Wface;

        Wface.rho =
            0.5 * (WL.rho + WR.rho);

        Wface.u =
            0.5 * (WL.u + WR.u);

        Wface.v =
            0.5 * (WL.v + WR.v);

        Wface.p =
            0.5 * (WL.p + WR.p);

        const Conserved Uface =
            primitive_to_conserved(Wface);

        return normal_flux(
            Uface,
            nx,
            ny
        );
    }

    void EulerSolver::validate_grid( const std::vector<std::vector<double>>& x_grid, const std::vector<std::vector<double>>& y_grid) const {
        if (x_grid.empty() || y_grid.empty())
            throw std::runtime_error("x_grid and y_grid must not be empty.");

        if (x_grid.size() != y_grid.size())
            throw std::runtime_error("x_grid and y_grid must have the same number of rows.");

        const size_t rows = x_grid.size();
        const size_t cols = x_grid[0].size();

        if (rows < 2 || cols < 2)
            throw std::runtime_error("Grid must contain at least 2x2 nodes.");

        for (size_t j = 0; j < rows; ++j)
        {
            if (x_grid[j].size() != cols || y_grid[j].size() != cols)
                throw std::runtime_error("All grid rows must have the same number of columns.");

            if (x_grid[j].size() != y_grid[j].size())
                throw std::runtime_error("x_grid and y_grid row sizes do not match.");
        }
    }

    const CellGeom& EulerSolver::cell_geom(int i, int j) const {
        return geom_[idx(i, j)];
    }

    double EulerSolver::signed_quad_area(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {
        return 0.5 * (
            x1 * y2 - y1 * x2 +
            x2 * y3 - y2 * x3 +
            x3 * y4 - y3 * x4 +
            x4 * y1 - y4 * x1
        );
    }

    FaceGeom EulerSolver::make_face( double xa, double ya, double xb, double yb, double cell_xc, double cell_yc) {
        FaceGeom f;

        const double dx = xb - xa;
        const double dy = yb - ya;

        f.length = std::sqrt(dx * dx + dy * dy);

        if (f.length <= 0.0)
            throw std::runtime_error("Zero-length face detected.");

        f.xc = 0.5 * (xa + xb);
        f.yc = 0.5 * (ya + yb);

        // Candidate outward normal for normal cell ordering.
        f.sx = dy;
        f.sy = -dx;

        // Check orientation using vector from cell center to face center.
        const double rx = f.xc - cell_xc;
        const double ry = f.yc - cell_yc;

        if (rx * f.sx + ry * f.sy < 0.0) {
            f.sx *= -1.0;
            f.sy *= -1.0;
        }

        f.nx = f.sx / f.length;
        f.ny = f.sy / f.length;

        return f;
    }

    void EulerSolver::build_geometry(const std::vector<std::vector<double>>& x_grid, const std::vector<std::vector<double>>& y_grid) {
        for (int j = j_start_; j < j_end_; ++j) {
            const int jm = j - j_start_;

            for (int i = i_start_; i < i_end_; ++i) {
                const int im = i - i_start_;

                const double x1 = x_grid[jm    ][im    ];
                const double y1 = y_grid[jm    ][im    ];

                const double x2 = x_grid[jm    ][im + 1];
                const double y2 = y_grid[jm    ][im + 1];

                const double x3 = x_grid[jm + 1][im + 1];
                const double y3 = y_grid[jm + 1][im + 1];

                const double x4 = x_grid[jm + 1][im    ];
                const double y4 = y_grid[jm + 1][im    ];

                CellGeom c;

                const double a_signed = signed_quad_area(
                    x1, y1,
                    x2, y2,
                    x3, y3,
                    x4, y4
                );

                c.area = std::abs(a_signed);

                if (c.area <= 0.0)
                    throw std::runtime_error("Cell with non-positive area detected.");

                c.xc = 0.25 * (x1 + x2 + x3 + x4);
                c.yc = 0.25 * (y1 + y2 + y3 + y4);

                c.face[0] = make_face(x1, y1, x2, y2, c.xc, c.yc);
                c.face[1] = make_face(x2, y2, x3, y3, c.xc, c.yc);
                c.face[2] = make_face(x3, y3, x4, y4, c.xc, c.yc);
                c.face[3] = make_face(x4, y4, x1, y1, c.xc, c.yc);

                geom_[idx(i, j)] = c;
            }
        }
    }

    void EulerSolver::check_geometry() const {
        double min_area = 1e300;
        double max_area = 0.0;
        double min_face = 1e300;
        double max_face = 0.0;

        for (int j = j_start_; j < j_end_; ++j) {

            for (int i = i_start_; i < i_end_; ++i) {

                const CellGeom& c = cell_geom(i, j);

                min_area = std::min(min_area, c.area);
                max_area = std::max(max_area, c.area);

                for (int f = 0; f < 4; ++f)
                {
                    const FaceGeom& face = c.face[f];

                    if (face.length <= 0.0)
                        throw std::runtime_error("Invalid face length.");

                    const double nmag = std::sqrt(face.nx * face.nx + face.ny * face.ny);

                    if (std::abs(nmag - 1.0) > 1e-12)
                        throw std::runtime_error("Face normal is not unit length.");

                    min_face = std::min(min_face, face.length);
                    max_face = std::max(max_face, face.length);
                }
            }
        }

        std::cout << "Geometry check passed:\n";
        std::cout << "  ni_nodes = " << ni_nodes_ << "\n";
        std::cout << "  nj_nodes = " << nj_nodes_ << "\n";
        std::cout << "  ni_cells = " << ni_cells_ << "\n";
        std::cout << "  nj_cells = " << nj_cells_ << "\n";
        std::cout << "  ni_total = " << ni_total_ << "\n";
        std::cout << "  nj_total = " << nj_total_ << "\n";
        std::cout << "  min_area = " << std::scientific << min_area << "\n";
        std::cout << "  max_area = " << std::scientific << max_area << "\n";
        std::cout << "  min_face = " << std::scientific << min_face << "\n";
        std::cout << "  max_face = " << std::scientific << max_face << "\n";
    }

    MeshGrid load_indexed_mesh(const std::string& filename) {
        struct Entry {
            int i;
            int j;
            double x;
            double y;
        };

        std::ifstream file(filename);

        if (!file) {
            throw std::runtime_error("Could not open mesh file: " + filename);
        }

        std::vector<Entry> entries;

        int i, j;
        double x, y;

        int imax = -1;
        int jmax = -1;

        while (file >> i >> j >> x >> y) {
            entries.push_back({i, j, x, y});
            imax = std::max(imax, i);
            jmax = std::max(jmax, j);
        }

        if (entries.empty()) {
            throw std::runtime_error("Mesh file is empty.");
        }

        const int ni_nodes = imax + 1;
        const int nj_nodes = jmax + 1;

        MeshGrid grid;

        grid.x.assign(nj_nodes, std::vector<double>(ni_nodes, 0.0));
        grid.y.assign(nj_nodes, std::vector<double>(ni_nodes, 0.0));

        std::vector<std::vector<bool>> seen(nj_nodes, std::vector<bool>(ni_nodes, false));

        for (const auto& e : entries) {

            if (e.i < 0 || e.i >= ni_nodes || e.j < 0 || e.j >= nj_nodes) {
                throw std::runtime_error("Invalid mesh index.");
            }

            if (seen[e.j][e.i]) {
                throw std::runtime_error("Duplicate mesh node.");
            }

            grid.x[e.j][e.i] = e.x;
            grid.y[e.j][e.i] = e.y;

            seen[e.j][e.i] = true;
        }

        for (int jj = 0; jj < nj_nodes; ++jj) {

            for (int ii = 0; ii < ni_nodes; ++ii) {

                if (!seen[jj][ii]) {
                    throw std::runtime_error("Missing node in structured mesh.");
                }
            }
        }

        std::cout << "Loaded indexed mesh:\n";
        std::cout << "  ni_nodes = " << ni_nodes << "\n";
        std::cout << "  nj_nodes = " << nj_nodes << "\n";

        return grid;
    }

    Conserved central_flux(const Conserved& UL, const Conserved& UR, double nx, double ny) {
        const Conserved FL = normal_flux(UL, nx, ny);

        const Conserved FR = normal_flux(UR, nx, ny);

        return 0.5 * (FL + FR);
    }

    double face_spectral_radius(const Conserved& UL, const Conserved& UR, double nx, double ny) {
        const Primitive WL = conserved_to_primitive(UL);

        const Primitive WR = conserved_to_primitive(UR);

        const double unL = WL.u * nx + WL.v * ny;

        const double unR = WR.u * nx + WR.v * ny;

        const double aL = sound_speed(WL);

        const double aR = sound_speed(WR);

        return 0.5 * (std::abs(unL) + aL + std::abs(unR) + aR);
    }

    double pressure_sensor(double p_minus, double p_center, double p_plus) {

        if (p_minus <= 0.0 || p_center <= 0.0 || p_plus <= 0.0) {

            throw std::runtime_error("Non-positive pressure in pressure_sensor.");
        }

        const double numerator = std::abs(p_plus - 2.0 * p_center + p_minus);

        const double denominator = p_plus + 2.0 * p_center + p_minus;

        return numerator / denominator;
    }

    Conserved jst_artificial_flux(const Conserved& U_minus, const Conserved& U_left, const Conserved& U_right, const Conserved& U_plus, double lambda_face, double epsilon_2, double epsilon_4) {

        Conserved D;

        D.rho = lambda_face * (epsilon_2 * (U_right.rho - U_left.rho) - epsilon_4 * (U_plus.rho - 3.0 * U_right.rho + 3.0 * U_left.rho - U_minus.rho ) );

        D.rhou = lambda_face * (epsilon_2 * (U_right.rhou - U_left.rhou) - epsilon_4 * (U_plus.rhou - 3.0 * U_right.rhou + 3.0 * U_left.rhou - U_minus.rhou ) );

        D.rhov = lambda_face * (epsilon_2 * (U_right.rhov - U_left.rhov) - epsilon_4 * (U_plus.rhov - 3.0 * U_right.rhov + 3.0 * U_left.rhov - U_minus.rhov ) );

        D.rhoE = lambda_face * (epsilon_2 * (U_right.rhoE - U_left.rhoE) - epsilon_4 * (U_plus.rhoE - 3.0 * U_right.rhoE + 3.0 * U_left.rhoE - U_minus.rhoE ) );

        return D;
    }

    static double pressure_from_state(const Conserved& U) {

        const Primitive W = conserved_to_primitive(U);

        return W.p;
    }

    void EulerSolver::compute_interior_periodic_residual_jst(const JSTParameters& parameters, bool recompute_dissipation) {

        zero_residual();

        if (recompute_dissipation) {
            for (auto& d : D_) {
                d = Conserved{};
            }
        }

        update_periodic_ghost_cells();

        // ------------------------------------------------------------
        // i-direction: periodic direction
        // ------------------------------------------------------------
        for (int j = j_start_; j < j_end_; ++j) {
            
            for (int i = i_start_; i < i_end_; ++i) {

                const int iL = i;

                int iR_physical = i + 1;

                if (iR_physical == i_end_) {

                    iR_physical = i_start_;
                }

                const CellGeom& cL = cell_geom(iL, j);

                const CellGeom& cR = cell_geom(iR_physical, j);

                const FaceGeom& face = cL.face[1];

                const Conserved& U_left = state(i, j);

                const Conserved& U_right = state(i + 1, j);

                // Central (convective) flux is recomputed on every RK stage.
                const Conserved Fc = central_flux(U_left, U_right, face.nx, face.ny);

                residual(iL, j) = residual(iL, j) + Fc * (face.length / cL.area);

                residual(iR_physical, j) = residual(iR_physical, j) - Fc * (face.length / cR.area);

                // Artificial dissipation is frozen across the RK stages.
                if (recompute_dissipation) {
                    const Conserved& U_minus = state(i - 1, j);

                    const Conserved& U_plus = state(i + 2, j);

                    double epsilon_2;
                    double epsilon_4;

                    if (parameters.first_order) {
                        // Coarse-grid first-order dissipation: constant scalar
                        // 2nd-difference everywhere, no 4th-difference. Skips
                        // the pressure sensor entirely (it can throw on a
                        // non-positive coarse-grid pressure, which is exactly
                        // the situation this robust operator exists to absorb).
                        epsilon_2 = 0.5;
                        epsilon_4 = 0.0;
                    }
                    else {
                        const double sensor_left = pressure_sensor( pressure_from_state(U_minus),pressure_from_state(U_left), pressure_from_state(U_right) );

                        const double sensor_right = pressure_sensor(pressure_from_state(U_left), pressure_from_state(U_right), pressure_from_state(U_plus));

                        epsilon_2 = parameters.k2 * std::max(sensor_left, sensor_right);

                        epsilon_4 = std::max(0.0, parameters.k4 - epsilon_2);
                    }

                    const double lambda_face = face_spectral_radius(U_left, U_right, face.nx, face.ny);

                    const Conserved D = jst_artificial_flux(U_minus, U_left, U_right, U_plus, lambda_face, epsilon_2, epsilon_4);

                    // F = Fc - D, so the dissipative part is -D on iL, +D on iR.
                    D_[idx(iL, j)] = D_[idx(iL, j)] - D * (face.length / cL.area);

                    D_[idx(iR_physical, j)] = D_[idx(iR_physical, j)] + D * (face.length / cR.area);
                }
            }
        }

        // ------------------------------------------------------------
        // j-direction: radial direction
        // ------------------------------------------------------------
        for (int j = j_start_; j < j_end_ - 1; ++j) {

            for (int i = i_start_; i < i_end_; ++i) {
                
                const int jB = j;
                const int jT = j + 1;

                const CellGeom& cB = cell_geom(i, jB);

                const CellGeom& cT = cell_geom(i, jT);

                const FaceGeom& face = cB.face[2];

                const Conserved& U_bottom = state(i, jB);

                const Conserved& U_top = state(i, jT);

                // Central (convective) flux is recomputed on every RK stage.
                const Conserved Fc = central_flux(U_bottom, U_top, face.nx, face.ny);

                residual(i, jB) = residual(i, jB) + Fc * (face.length / cB.area);

                residual(i, jT) = residual(i, jT) - Fc * (face.length / cT.area);

                if (!recompute_dissipation) {
                    continue;
                }

                const double lambda_face = face_spectral_radius(U_bottom, U_top, face.nx, face.ny);

                Conserved U_bottom_bottom;
                Conserved U_top_top;

                // --------------------------------------------------------
                // Build the missing JST stencil states near boundaries.
                //
                // Interior face:
                //     U_bottom_bottom = U(j-1)
                //     U_top_top       = U(j+2)
                //
                // Wall-adjacent face:
                //     U_bottom_bottom = reflected wall ghost
                //
                // Farfield-adjacent face:
                //     U_top_top = zero-gradient ghost = U_top
                // --------------------------------------------------------

                if (j - 1 >= j_start_) {
                    U_bottom_bottom = state(i, j - 1);
                }

                else {
                    const CellGeom& c_wall = cell_geom(i, jB);

                    const FaceGeom& wall_face = c_wall.face[0];

                    const Primitive W_bottom = conserved_to_primitive(U_bottom);

                    const Primitive W_ghost = reflected_wall_state(W_bottom, wall_face.nx, wall_face.ny);

                    U_bottom_bottom = primitive_to_conserved(W_ghost);
                }

                if (j + 2 < j_end_) {
                    U_top_top = state(i, j + 2);
                }

                else {
                    // Farfield-adjacent JST stencil closure.
                    // Simple zero-gradient extrapolation.
                    U_top_top =  U_top;
                }

                double epsilon_2;
                double epsilon_4;

                if (parameters.first_order) {
                    // Coarse-grid first-order dissipation (see i-direction).
                    epsilon_2 = 0.5;
                    epsilon_4 = 0.0;
                }
                else {
                    const double sensor_bottom = pressure_sensor(pressure_from_state(U_bottom_bottom), pressure_from_state(U_bottom), pressure_from_state(U_top));

                    const double sensor_top = pressure_sensor(pressure_from_state(U_bottom), pressure_from_state(U_top), pressure_from_state(U_top_top));

                    epsilon_2 = parameters.k2 * std::max(sensor_bottom, sensor_top);

                    epsilon_4 = std::max(0.0, parameters.k4 - epsilon_2);
                }

                const Conserved D = jst_artificial_flux(U_bottom_bottom, U_bottom, U_top, U_top_top, lambda_face, epsilon_2, epsilon_4);

                // F = Fc - D, so the dissipative part is -D on jB, +D on jT.
                D_[idx(i, jB)] = D_[idx(i, jB)] - D * (face.length / cB.area);

                D_[idx(i, jT)] = D_[idx(i, jT)] + D * (face.length / cT.area);
            }
        }
    }

    void EulerSolver::compute_full_residual_jst(const Primitive& Winf, const JSTParameters& parameters, bool recompute_dissipation) {
        compute_interior_periodic_residual_jst(parameters, recompute_dissipation);


        add_wall_boundary_residual();

        add_farfield_boundary_residual(Winf);

        // Add the (possibly frozen) artificial dissipation to the convective
        // and boundary residual already accumulated in R_.
        for (int j = j_start_; j < j_end_; ++j) {
            for (int i = i_start_; i < i_end_; ++i) {
                const int id = idx(i, j);
                R_[id] = R_[id] + D_[id];
            }
        }
    }

    void EulerSolver::compute_irs_coefficients(double eps) {

        // Uniform coefficient in both directions. Smoothing the RK increment
        // dt*R (not R) means a large eps over-couples a grid line and
        // destabilizes; mild eps ~ 0.1 is stable to CFL ~ 8 on the stretched mesh.
        for (int j = j_start_; j < j_end_; ++j) {
            for (int i = i_start_; i < i_end_; ++i) {
                eps_i_[idx(i, j)] = eps;
                eps_j_[idx(i, j)] = eps;
            }
        }
    }

    void EulerSolver::smooth_residual_implicit(const std::vector<double>& dt) {

        const int ni = ni_cells_;
        const int nj = nj_cells_;
        const int nmax = std::max(ni, nj);

        std::vector<double> sub(nmax), diag(nmax), sup(nmax), rhs(nmax);
        std::vector<double> bb(nmax), uvec(nmax), cprime(nmax);

        // Per-component accessor for a Conserved residual.
        auto comp = [](Conserved& C, int c) -> double& {
            switch (c) {
                case 0:  return C.rho;
                case 1:  return C.rhou;
                case 2:  return C.rhov;
                default: return C.rhoE;
            }
        };

        // We smooth the update increment delta = dt * R (not R itself). On a
        // mesh with a ~3.5-order cell-size spread, every cell's increment is
        // O(1) regardless of size, so coupling a whole grid line cannot leak
        // an outsized update between dissimilar cells. We multiply by dt on
        // gather and divide back on scatter so the existing RK update
        // U = U0 - alpha*dt*R_ reproduces U0 - alpha*(smoothed increment).

        // ---- i-direction sweep: periodic around the airfoil ----
        for (int j = j_start_; j < j_end_; ++j) {

            for (int m = 0; m < ni; ++m) {
                const double e = eps_i_[idx(i_start_ + m, j)];
                sub[m]  = -e;
                sup[m]  = -e;
                diag[m] = 1.0 + 2.0 * e;
            }

            for (int c = 0; c < 4; ++c) {

                for (int m = 0; m < ni; ++m) {
                    const int id = idx(i_start_ + m, j);
                    rhs[m] = comp(R_[id], c) * dt[id];
                }

                cyclic_scalar(sub, diag, sup, rhs, ni, bb, uvec, cprime);

                for (int m = 0; m < ni; ++m) {
                    const int id = idx(i_start_ + m, j);
                    comp(R_[id], c) = rhs[m] / dt[id];
                }
            }
        }

        // ---- j-direction sweep: Neumann at wall (j_start) and farfield ----
        for (int i = i_start_; i < i_end_; ++i) {

            for (int m = 0; m < nj; ++m) {
                const double e = eps_j_[idx(i, j_start_ + m)];
                const bool first = (m == 0);
                const bool last  = (m == nj - 1);

                sub[m] = first ? 0.0 : -e;
                sup[m] = last  ? 0.0 : -e;

                // Zero-gradient ghost: a boundary cell keeps only its single
                // physical neighbor, so its diagonal is 1 + e instead of 1 + 2e.
                const int n_off = (first ? 0 : 1) + (last ? 0 : 1);
                diag[m] = 1.0 + n_off * e;
            }

            for (int c = 0; c < 4; ++c) {

                for (int m = 0; m < nj; ++m) {
                    const int id = idx(i, j_start_ + m);
                    rhs[m] = comp(R_[id], c) * dt[id];
                }

                thomas_scalar(sub, diag, sup, rhs, nj, cprime);

                for (int m = 0; m < nj; ++m) {
                    const int id = idx(i, j_start_ + m);
                    comp(R_[id], c) = rhs[m] / dt[id];
                }
            }
        }
    }

    void EulerSolver::set_uniform_state(const Primitive& W) {
        
        const Conserved U = primitive_to_conserved(W); 

        for (int j = j_start_; j < j_end_; ++j) { 

            for (int i = i_start_; i < i_end_; ++i) { 

                U_[idx(i, j)] = U;
            }
        }
    }

    const Conserved& EulerSolver::state(int i, int j) const { 
        return U_[idx(i, j)];
    }

    Conserved& EulerSolver::state(int i, int j) { 
        return U_[idx(i, j)]; 
    }

    void EulerSolver::update_periodic_ghost_cells() { 

        for (int j = j_start_; j < j_end_; ++j) { 

            for (int g = 0 ; g < ng_; ++g) { 

                const int i_left_ghost = i_start_ - ng_ + g; 
                const int i_left_src = i_end_ - ng_ + g; 

                const int i_right_ghost = i_end_ + g; 
                const int i_right_src = i_start_ + g; 

                U_[idx(i_left_ghost, j)] = U_[idx(i_left_src, j)]; 
                U_[idx(i_right_ghost, j)] = U_[idx(i_right_src, j)];

            }
        }

    }

    const Conserved& EulerSolver::residual(int i, int j) const {
        return R_[idx(i, j)];
    }

    Conserved& EulerSolver::residual(int i, int j) {
        return R_[idx(i, j)];
    }

    void EulerSolver::zero_residual() {
        for (auto& r : R_)
        {
            r = Conserved{};
        }
    }


    Primitive reflected_wall_state(const Primitive& Wi, double nx, double ny) { 

        if (!is_physical(Wi)) { 

            throw std::runtime_error("Non-physical state in reflected_wall_state.");

        }

        const double un = Wi.u * nx + Wi.v * ny; 

        Primitive Wg; 
        Wg.rho = Wi.rho; 
        Wg.p = Wi.p; 

        Wg.u = Wi.u - 2.0 * un * nx; 
        Wg.v = Wi.v - 2.0 * un * ny; 

        return Wg; 

    }

    Primitive farfield_riemann_state(const Primitive& Wi,const Primitive& Winf,double nx,double ny) {

        if (!is_physical(Wi)) {
            throw std::runtime_error("Non-physical interior state in farfield_riemann_state.");
        }

        if (!is_physical(Winf)) {
            throw std::runtime_error("Non-physical freestream state in farfield_riemann_state.");
        }

        const double rhoi = Wi.rho;
        const double ui   = Wi.u;
        const double vi   = Wi.v;
        const double pi   = Wi.p;

        const double rhoinf = Winf.rho;
        const double uinf   = Winf.u;
        const double vinf   = Winf.v;
        const double pinf   = Winf.p;

        const double ai = sound_speed(Wi);
        const double ainf = sound_speed(Winf);

        const double uni = ui * nx + vi * ny;
        const double uninf = uinf * nx + vinf * ny;

        const double Vi = std::sqrt(ui * ui + vi * vi);
        const double Mach_i = Vi / ai;

        // Supersonic cases.
        // Normal points outward from the computational domain.
        // uni < 0 means flow enters the domain.
        if (Mach_i >= 1.0) {

            if (uni < 0.0) {
                // Supersonic inflow: all information enters from freestream.
                return Winf;
            }
            else {
                // Supersonic outflow: all information comes from interior.
                return Wi;
            }
        }

        // Tangential velocity convention:
        // t = (-ny, nx)
        // u_t = velocity dot t = -u*ny + v*nx
        if (uni < 0.0) {
            // Subsonic inflow:
            // Incoming information mostly from freestream,
            // outgoing acoustic information from interior.
            const double pb = 0.5 * (pinf + pi - rhoi * ai * (uninf - uni));

            const double rhob = rhoinf + (pb - pinf) / (ainf * ainf);

            const double unb = uninf - (pinf - pb) / (rhoinf * ainf);

            const double ut_inf = -uinf * ny + vinf * nx;

            Primitive Wb;
            Wb.rho = rhob;
            Wb.p   = pb;

            Wb.u = unb * nx - ut_inf * ny;
            Wb.v = unb * ny + ut_inf * nx;

            return Wb;
        }

        else {
            // Subsonic outflow:
            // Pressure is imposed from farfield,
            // velocity/density mostly extrapolated from interior.
            const double pb = pinf;

            const double rhob = rhoi + (pb - pi) / (ai * ai);

            const double unb = uni + (pi - pb) / (rhoi * ai);

            const double ut_i = -ui * ny + vi * nx;

            Primitive Wb;
            Wb.rho = rhob;
            Wb.p   = pb;

            Wb.u = unb * nx - ut_i * ny;
            Wb.v = unb * ny + ut_i * nx;

            return Wb;
        }
    }

    void EulerSolver::add_wall_boundary_residual() {

        const int j = j_start_;

        for (int i = i_start_; i < i_end_; ++i) {

            const CellGeom& c = cell_geom(i, j);

            const FaceGeom& face = c.face[0];

            const Conserved& Ui = state(i, j);

            const Primitive Wi = conserved_to_primitive(Ui);

            const Primitive Wg = reflected_wall_state(Wi, face.nx, face.ny);

            Primitive Wface;

            Wface.rho = 0.5 * (Wi.rho + Wg.rho);

            Wface.u = 0.5 * (Wi.u + Wg.u);

            Wface.v = 0.5 * (Wi.v + Wg.v);

            Wface.p = 0.5 * (Wi.p + Wg.p);

            const Conserved Uface = primitive_to_conserved(Wface);

            const Conserved F = normal_flux(Uface, face.nx, face.ny);

            residual(i, j) = residual(i, j) + F * (face.length / c.area);
        }
    }

    void EulerSolver::add_farfield_boundary_residual(const Primitive& Winf) {

        const int j = j_end_ - 1;

        for (int i = i_start_; i < i_end_; ++i) {

            const CellGeom& c = cell_geom(i, j);

            const FaceGeom& face = c.face[2]; // top farfield face

            const Conserved& Ui = state(i, j);

            const Primitive Wi = conserved_to_primitive(Ui);

            const Primitive Wb = farfield_riemann_state(Wi, Winf, face.nx, face.ny);

            const Conserved F = averaged_state_flux(Wi, Wb, face.nx, face.ny); 

            residual(i, j) = residual(i, j) + F * (face.length / c.area);
        }
    }

    double EulerSolver::local_time_step(int i, int j, double cfl) const { 

        if (cfl <= 0) { 

            throw std::runtime_error("CFL must be positive"); 

        }

        const CellGeom& c = cell_geom(i, j); 
        const Primitive W = conserved_to_primitive(state(i,j)); 

        const double a = sound_speed(W); 

        double spectral_sum = 0.0;

        for (int f = 0; f < 4; ++f) { 

            const FaceGeom face = c.face[f]; 

            const double un = W.u * face.nx + W.v * face.ny; 

            spectral_sum += (std::abs(un) + a) * face.length; 

        }

        if (spectral_sum <= 0.0) { 
            throw std::runtime_error("Non-positive spectral sum in local_time_step.");
        }

        return cfl * c.area / spectral_sum; 

    }

    std::vector<double> EulerSolver::compute_local_time_steps(double cfl) const { 

        std::vector<double> dt(ni_total_* nj_total_, 0.0); 

        for (int j = j_start_; j < j_end_; ++j) { 

            for (int i = i_start_; i < i_end_; ++i) { 

                dt[idx(i, j)] = local_time_step(i, j, cfl); 

            }

        }

        return dt; 

    }

    double EulerSolver::residual_linf_current() const { 

        double max_value = 0.0; 

        for (int j = j_start_; j < j_end_; ++j) { 

            for (int i = i_start_; i < i_end_; ++i) { 

                const Conserved& r = residual(i, j); 

                max_value = std::max(max_value, std::abs(r.rho)); 
                max_value = std::max(max_value, std::abs(r.rhou)); 
                max_value = std::max(max_value, std::abs(r.rhov)); 
                max_value = std::max(max_value, std::abs(r.rhoE)); 

            }
        }

        return max_value; 

    }

    Primitive freestream_state(double mach, double alpha_rad) {
        Primitive W;

        W.rho = 1.0;
        W.p   = 1.0 / GAMMA;

        const double a_inf = sound_speed(W);

        const double V_inf = mach * a_inf;

        W.u =  V_inf * std::cos(alpha_rad);

        // Positive alpha convention:
        // flow comes from below the chord line in the body-fixed frame,
        // giving positive lift for positive alpha (matches Vassberg-Jameson).
        W.v =  V_inf * std::sin(alpha_rad);

        return W;
    }

    void EulerSolver::rk5_pseudo_step(const Primitive& Winf, const JSTParameters& parameters, double cfl_now) {

        static const std::array<double, 5> alpha =
        {
            1.0 / 4.0,
            1.0 / 6.0,
            3.0 / 8.0,
            1.0 / 2.0,
            1.0
        };

        // Save U^n. Each RK stage is built from this reference state.
        const std::vector<Conserved> U0 = U_;

        // Local time steps for this iteration, frozen across the RK stages.
        const std::vector<double> dt = compute_local_time_steps(cfl_now);

        // Implicit-residual-smoothing coefficient, frozen across the RK stages
        // like dt. Mild uniform value; stable to CFL ~ 8.
        const double irs_eps = 0.1;
        compute_irs_coefficients(irs_eps);

        for (int stage = 0; stage < 5; ++stage) {
            // Convective + boundary fluxes are recomputed every stage, but the
            // JST artificial dissipation is frozen: the caller evaluated the
            // residual with recompute_dissipation=true before this call, so D_
            // already holds the dissipation and each stage reuses it.
            compute_full_residual_jst(Winf, parameters, /*recompute_dissipation=*/false);

            // FAS forcing: relax toward R = tau by driving the modified
            // residual R - forcing_ to zero. forcing_ is zero on the fine
            // (single) grid, so this is an exact no-op there (x - 0 == x) and
            // the single-grid trajectory is unchanged.
            for (int j = j_start_; j < j_end_; ++j) {
                for (int i = i_start_; i < i_end_; ++i) {
                    const int id = idx(i, j);
                    R_[id] = R_[id] - forcing_[id];
                }
            }

            // Implicit residual smoothing extends the stable CFL limit of the
            // 5-stage scheme. Applied only to the RK update residual; the
            // residual reported by the driver remains unsmoothed.
            smooth_residual_implicit(dt);

            const double a_stage = alpha[stage];

            for (int j = j_start_; j < j_end_; ++j) {

                for (int i = i_start_; i < i_end_; ++i) {

                    const int id = idx(i, j);

                    U_[id] =
                        U0[id] - a_stage * dt[id] * R_[id];

                    const Primitive Wnew = conserved_to_primitive(U_[id]);

                    if (!is_physical(Wnew)) {
                        throw std::runtime_error("Non-physical state during RK5 pseudo-time iteration.");
                    }
                }
            }

            update_periodic_ghost_cells();
        }
    }

    void EulerSolver::write_flowfield(const std::string& filename) const {

        std::ofstream file(filename);

        if (!file) {
            throw std::runtime_error("Could not open flowfield output file: " + filename);
        }

        file << "# i j xc yc rho u v p mach\n";
        file << std::scientific << std::setprecision(12);

        for (int j = j_start_; j < j_end_; ++j) {

            for (int i = i_start_; i < i_end_; ++i) {

                const CellGeom& c = cell_geom(i, j);
                const Conserved& U = state(i, j);
                const Primitive W = conserved_to_primitive(U);

                const double a = sound_speed(W);
                const double V = std::sqrt(W.u * W.u + W.v * W.v);
                const double mach = V / a;

                // Output physical-cell indices, not ghost-offset indices.
                const int im = i - i_start_;
                const int jm = j - j_start_;

                file
                    << im << " "
                    << jm << " "
                    << c.xc << " "
                    << c.yc << " "
                    << W.rho << " "
                    << W.u << " "
                    << W.v << " "
                    << W.p << " "
                    << mach << "\n";
            }
        }

        std::cout << "Wrote flowfield: " << filename << "\n";
    }

    void EulerSolver::write_surface_cp(const std::string& filename, const Primitive& Winf) const {

        std::ofstream out(filename);

        if (!out) {
            throw std::runtime_error("Could not open surface Cp file: " + filename);
        }

        const double Vinf2 = Winf.u * Winf.u + Winf.v * Winf.v;

        const double qinf = 0.5 * Winf.rho * Vinf2;

        if (qinf <= 0.0) {
            throw std::runtime_error("Invalid dynamic pressure in write_surface_cp.");
        }

        out << "# i x_wall y_wall p_first cp_first p_wall_extrapolated cp_wall_extrapolated\n";
        out << std::scientific << std::setprecision(12);

        const int j0 = j_start_;

        const int j1 = j_start_ + 1;

        if (j1 >= j_end_) {
            throw std::runtime_error("Not enough radial cells for wall pressure extrapolation.");
        }

        for (int i = i_start_; i < i_end_; ++i) {
            const CellGeom& c = cell_geom(i, j0);

            
            const FaceGeom& wall_face = c.face[0];

            const Primitive W0 = conserved_to_primitive(state(i, j0));

            const Primitive W1 = conserved_to_primitive(state(i, j1));

            // First-cell estimate.
            const double p_first = W0.p;

            const double cp_first = (p_first - Winf.p) / qinf;

            // Linear extrapolation from the first two cell centers toward the wall:
            //
            //     p_wall ≈ 1.5 p_0 - 0.5 p_1
            //
            // This assumes the first two cell centers are approximately equally
            // spaced in the wall-normal direction.
            const double p_wall = 1.5 * W0.p - 0.5 * W1.p;

            const double cp_wall = (p_wall - Winf.p) / qinf;

            out << (i - i_start_) << " "
                << wall_face.xc << " "
                << wall_face.yc << " "
                << p_first << " "
                << cp_first << " "
                << p_wall << " "
                << cp_wall << "\n";
        }

        std::cout << "Wrote surface Cp: "
                << filename
                << "\n";
    }

    AeroCoefficients EulerSolver::compute_aero_coefficients(
        const Primitive& Winf, double alpha_rad) const
    {
        const double Vinf2 = Winf.u * Winf.u + Winf.v * Winf.v;
        const double qinf  = 0.5 * Winf.rho * Vinf2;

        if (qinf <= 0.0) {
            throw std::runtime_error(
                "Invalid dynamic pressure in compute_aero_coefficients.");
        }

        constexpr double chord_ref = 1.0;
        constexpr double x_ref = 0.25; // quarter-chord moment reference
        constexpr double y_ref = 0.0;

        const double cosa = std::cos(alpha_rad);
        const double sina = std::sin(alpha_rad);

        double cx = 0.0, cy = 0.0, cm = 0.0;
        const int j = j_start_; // wall row

        for (int i = i_start_; i < i_end_; ++i) {
            const FaceGeom& face = cell_geom(i, j).face[0]; // wall face
            const Primitive W = conserved_to_primitive(state(i, j));
            const double cp = (W.p - Winf.p) / qinf;

            // face.n is the fluid cell's outward normal (into the body); the
            // body normal is -face.n and pressure acts along -body_normal, so
            // the force on the body is +cp*face.n (the two signs cancel).
            const double dCx = cp * face.nx * face.length / chord_ref;
            const double dCy = cp * face.ny * face.length / chord_ref;
            cx += dCx;
            cy += dCy;

            // Moment about (x_ref, y_ref): dCm = (x_ref-x)*dCy - (y_ref-y)*dCx.
            cm += (x_ref - face.xc) * dCy - (y_ref - face.yc) * dCx;
        }

        const double cd = cx * cosa + cy * sina;
        const double cl = -cx * sina + cy * cosa;

        return AeroCoefficients{
            cl,
            cd,
            cm
        };
    }

}