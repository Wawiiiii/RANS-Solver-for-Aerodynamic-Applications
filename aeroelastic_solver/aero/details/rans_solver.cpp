#include "rans_solver.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace rans {

    Vec2 least_squares_gradient(
        const double dx[], const double dy[], const double dphi[], int n)
    {
        double a11 = 0.0, a12 = 0.0, a22 = 0.0, b1 = 0.0, b2 = 0.0;

        for (int k = 0; k < n; ++k) {
            const double w = 1.0 / (dx[k] * dx[k] + dy[k] * dy[k]);
            a11 += w * dx[k] * dx[k];
            a12 += w * dx[k] * dy[k];
            a22 += w * dy[k] * dy[k];
            b1  += w * dx[k] * dphi[k];
            b2  += w * dy[k] * dphi[k];
        }

        const double det = a11 * a22 - a12 * a12;
        if (det == 0.0)
            throw std::runtime_error("rans: degenerate least-squares stencil.");

        return Vec2{
            ( a22 * b1 - a12 * b2) / det,
            (-a12 * b1 + a11 * b2) / det
        };
    }

    Vec2 corrected_face_gradient(
        const Vec2& grad_L, const Vec2& grad_R,
        double phi_L, double phi_R,
        double dx_LR, double dy_LR)
    {
        const double d = std::sqrt(dx_LR * dx_LR + dy_LR * dy_LR);
        const double ex = dx_LR / d;
        const double ey = dy_LR / d;

        const Vec2 avg{
            0.5 * (grad_L.x + grad_R.x),
            0.5 * (grad_L.y + grad_R.y)
        };

        const double avg_along_e = avg.x * ex + avg.y * ey;
        const double dphi_de = (phi_R - phi_L) / d;
        const double correction = dphi_de - avg_along_e;

        return Vec2{
            avg.x + correction * ex,
            avg.y + correction * ey
        };
    }

    Conserved viscous_normal_flux(
        const Vec2& grad_u, const Vec2& grad_v, const Vec2& grad_T,
        double u_face, double v_face,
        double mu, double k,
        double nx, double ny)
    {
        const double du_dx = grad_u.x;
        const double du_dy = grad_u.y;
        const double dv_dx = grad_v.x;
        const double dv_dy = grad_v.y;

        const double divergence = du_dx + dv_dy;

        const double tau_xx = mu * (2.0 * du_dx - (2.0 / 3.0) * divergence);
        const double tau_yy = mu * (2.0 * dv_dy - (2.0 / 3.0) * divergence);
        const double tau_xy = mu * (du_dy + dv_dx);

        Conserved Fv;

        Fv.rho = 0.0;
        Fv.rhou = tau_xx * nx + tau_xy * ny;
        Fv.rhov = tau_xy * nx + tau_yy * ny;

        // energy: work by stresses (V.tau) plus conduction (+k grad T, since q = -k grad T)
        const double theta_x = u_face * tau_xx + v_face * tau_xy + k * grad_T.x;
        const double theta_y = u_face * tau_xy + v_face * tau_yy + k * grad_T.y;

        Fv.rhoE = theta_x * nx + theta_y * ny;

        return Fv;
    }

    // ------------------------------------------------------------------
    // State conversions.
    // ------------------------------------------------------------------

    bool is_physical(const Primitive& W)
    {
        return W.rho > 0.0 && W.p > 0.0;
    }

    Primitive conserved_to_primitive(const Conserved& U)
    {
        if (U.rho <= 0.0)
            throw std::runtime_error("rans: non-positive density.");

        Primitive W;
        W.rho = U.rho;
        W.u = U.rhou / U.rho;
        W.v = U.rhov / U.rho;

        const double E = U.rhoE / U.rho;
        W.p = (GAMMA - 1.0) * U.rho * (E - 0.5 * (W.u * W.u + W.v * W.v));

        if (W.p <= 0.0)
            throw std::runtime_error("rans: non-positive pressure.");

        return W;
    }

    Conserved primitive_to_conserved(const Primitive& W)
    {
        if (!is_physical(W))
            throw std::runtime_error("rans: non-physical primitive state.");

        const double E =
            W.p / ((GAMMA - 1.0) * W.rho) + 0.5 * (W.u * W.u + W.v * W.v);

        Conserved U;
        U.rho = W.rho;
        U.rhou = W.rho * W.u;
        U.rhov = W.rho * W.v;
        U.rhoE = W.rho * E;
        return U;
    }

    double sound_speed(const Primitive& W)
    {
        if (!is_physical(W))
            throw std::runtime_error("rans: non-physical primitive state in sound_speed.");

        return std::sqrt(GAMMA * W.p / W.rho);
    }

    double temperature(const Primitive& W)
    {
        return W.p / W.rho;
    }

    // ------------------------------------------------------------------
    // Convective / inviscid mean-flow flux atoms.
    // ------------------------------------------------------------------

    Conserved normal_flux(const Conserved& U, double nx, double ny)
    {
        const Primitive W = conserved_to_primitive(U);

        const double un = W.u * nx + W.v * ny;

        Conserved F;
        F.rho  = W.rho * un;
        F.rhou = W.rho * W.u * un + W.p * nx;
        F.rhov = W.rho * W.v * un + W.p * ny;
        F.rhoE = (U.rhoE + W.p) * un;
        return F;
    }

    Conserved central_flux(const Conserved& UL, const Conserved& UR, double nx, double ny)
    {
        const Conserved FL = normal_flux(UL, nx, ny);
        const Conserved FR = normal_flux(UR, nx, ny);

        Conserved F;
        F.rho  = 0.5 * (FL.rho  + FR.rho);
        F.rhou = 0.5 * (FL.rhou + FR.rhou);
        F.rhov = 0.5 * (FL.rhov + FR.rhov);
        F.rhoE = 0.5 * (FL.rhoE + FR.rhoE);
        return F;
    }

    double face_spectral_radius(const Conserved& UL, const Conserved& UR, double nx, double ny)
    {
        const Primitive WL = conserved_to_primitive(UL);
        const Primitive WR = conserved_to_primitive(UR);

        const double unL = WL.u * nx + WL.v * ny;
        const double unR = WR.u * nx + WR.v * ny;

        return 0.5 * (
            std::abs(unL) + sound_speed(WL) +
            std::abs(unR) + sound_speed(WR));
    }

    // ------------------------------------------------------------------
    // Phase B0: mesh handling.
    // ------------------------------------------------------------------

    RansSolver::RansSolver(
        const std::vector<double>& x_nodes,
        const std::vector<double>& y_nodes,
        int n_t_nodes,
        int n_r_nodes)
    {
        ni_nodes_ = n_t_nodes;
        nj_nodes_ = n_r_nodes;

        validate_grid(x_nodes, y_nodes);

        ni_cells_ = ni_nodes_ - 1;
        nj_cells_ = nj_nodes_ - 1;

        ni_total_ = ni_cells_ + 2 * ng_;
        nj_total_ = nj_cells_ + 2 * ng_;

        i_start_ = ng_;
        i_end_   = ng_ + ni_cells_;
        j_start_ = ng_;
        j_end_   = ng_ + nj_cells_;

        const size_t total = static_cast<size_t>(ni_total_) * nj_total_;
        geom_.resize(total);
        U_.resize(total);
        grad_u_.resize(total);
        grad_v_.resize(total);
        grad_T_.resize(total);
        R_.resize(total);

        build_geometry(x_nodes, y_nodes);
    }

    void RansSolver::validate_grid(
        const std::vector<double>& x_nodes,
        const std::vector<double>& y_nodes) const
    {
        if (ni_nodes_ < 2 || nj_nodes_ < 2)
            throw std::runtime_error("RansSolver: need at least 2x2 nodes.");

        const size_t expected =
            static_cast<size_t>(ni_nodes_) * static_cast<size_t>(nj_nodes_);

        if (x_nodes.size() != expected || y_nodes.size() != expected)
            throw std::runtime_error("RansSolver: node array size does not match n_t*n_r.");
    }

    double RansSolver::signed_quad_area(
        double x1, double y1, double x2, double y2,
        double x3, double y3, double x4, double y4)
    {
        return 0.5 * (
            x1 * y2 - y1 * x2 +
            x2 * y3 - y2 * x3 +
            x3 * y4 - y3 * x4 +
            x4 * y1 - y4 * x1);
    }

    FaceGeom RansSolver::make_face(
        double xa, double ya, double xb, double yb,
        double cell_xc, double cell_yc)
    {
        FaceGeom f;

        const double dx = xb - xa;
        const double dy = yb - ya;

        f.length = std::sqrt(dx * dx + dy * dy);

        if (f.length <= 0.0)
            throw std::runtime_error("RansSolver: zero-length face.");

        f.xc = 0.5 * (xa + xb);
        f.yc = 0.5 * (ya + yb);

        // Rotate the edge to a candidate normal, then flip it outward if it
        // points back toward the cell center.
        f.sx = dy;
        f.sy = -dx;

        const double rx = f.xc - cell_xc;
        const double ry = f.yc - cell_yc;

        if (rx * f.sx + ry * f.sy < 0.0) {
            f.sx = -f.sx;
            f.sy = -f.sy;
        }

        f.nx = f.sx / f.length;
        f.ny = f.sy / f.length;

        return f;
    }

    void RansSolver::build_geometry(
        const std::vector<double>& x_nodes,
        const std::vector<double>& y_nodes)
    {
        auto node = [&](int im, int jm) {
            return static_cast<size_t>(im) + static_cast<size_t>(jm) * ni_nodes_;
        };

        for (int j = j_start_; j < j_end_; ++j) {
            const int jm = j - j_start_;

            for (int i = i_start_; i < i_end_; ++i) {
                const int im = i - i_start_;

                const double x1 = x_nodes[node(im,     jm)];
                const double y1 = y_nodes[node(im,     jm)];
                const double x2 = x_nodes[node(im + 1, jm)];
                const double y2 = y_nodes[node(im + 1, jm)];
                const double x3 = x_nodes[node(im + 1, jm + 1)];
                const double y3 = y_nodes[node(im + 1, jm + 1)];
                const double x4 = x_nodes[node(im,     jm + 1)];
                const double y4 = y_nodes[node(im,     jm + 1)];

                CellGeom c;

                const double a_signed =
                    signed_quad_area(x1, y1, x2, y2, x3, y3, x4, y4);

                c.area = std::abs(a_signed);

                if (c.area <= 0.0)
                    throw std::runtime_error("RansSolver: non-positive cell area.");

                c.xc = 0.25 * (x1 + x2 + x3 + x4);
                c.yc = 0.25 * (y1 + y2 + y3 + y4);

                c.face[0] = make_face(x1, y1, x2, y2, c.xc, c.yc); // j-
                c.face[1] = make_face(x2, y2, x3, y3, c.xc, c.yc); // i+
                c.face[2] = make_face(x3, y3, x4, y4, c.xc, c.yc); // j+
                c.face[3] = make_face(x4, y4, x1, y1, c.xc, c.yc); // i-

                geom_[idx(i, j)] = c;
            }
        }
    }

    GeometryReport RansSolver::check_geometry() const
    {
        GeometryReport r;
        r.ni_cells = ni_cells_;
        r.nj_cells = nj_cells_;
        r.min_area = 1e300;
        r.max_area = 0.0;
        r.min_face = 1e300;
        r.max_face = 0.0;

        for (int j = j_start_; j < j_end_; ++j) {
            for (int i = i_start_; i < i_end_; ++i) {
                const CellGeom& c = cell_geom(i, j);

                if (c.area <= 0.0) ++r.inverted_cells;

                r.min_area = std::min(r.min_area, c.area);
                r.max_area = std::max(r.max_area, c.area);

                for (int f = 0; f < 4; ++f) {
                    const FaceGeom& face = c.face[f];
                    r.min_face = std::min(r.min_face, face.length);
                    r.max_face = std::max(r.max_face, face.length);

                    const double nmag =
                        std::sqrt(face.nx * face.nx + face.ny * face.ny);
                    r.max_unit_normal_error =
                        std::max(r.max_unit_normal_error, std::abs(nmag - 1.0));
                }
            }
        }

        return r;
    }

    // ------------------------------------------------------------------
    // Phase B1: ghost filling + whole-field gradient pass.
    // ------------------------------------------------------------------

    void RansSolver::fill_ghost_cells()
    {
        // Periodic in i: copy the physical columns across the O-grid seam.
        for (int j = j_start_; j < j_end_; ++j)
            for (int g = 0; g < ng_; ++g) {
                U_[idx(i_start_ - ng_ + g, j)] = U_[idx(i_end_ - ng_ + g, j)];
                U_[idx(i_end_ + g, j)]         = U_[idx(i_start_ + g, j)];
            }

        // Zero-gradient in j at the wall and farfield, over the full padded
        // width so corner ghosts are consistent.
        for (int i = 0; i < ni_total_; ++i)
            for (int g = 1; g <= ng_; ++g) {
                U_[idx(i, j_start_ - g)]     = U_[idx(i, j_start_)];
                U_[idx(i, j_end_ - 1 + g)]   = U_[idx(i, j_end_ - 1)];
            }
    }

    void RansSolver::compute_gradients()
    {
        for (int j = j_start_; j < j_end_; ++j)
            for (int i = i_start_; i < i_end_; ++i) {
                const CellGeom& c = cell_geom(i, j);
                const Primitive Wc = conserved_to_primitive(U_[idx(i, j)]);
                const double Tc = temperature(Wc);

                double dx[4], dy[4], du[4], dv[4], dT[4];
                int n = 0;

                // Gather real cell neighbors. i wraps periodically across the
                // O-grid seam; j is one-sided at the wall / farfield (no ghost
                // geometry needed -- least squares stays exact with 3 neighbors).
                auto add = [&](int ni, int nj) {
                    const CellGeom& nc = cell_geom(ni, nj);
                    const Primitive Wn = conserved_to_primitive(U_[idx(ni, nj)]);
                    dx[n] = nc.xc - c.xc;
                    dy[n] = nc.yc - c.yc;
                    du[n] = Wn.u - Wc.u;
                    dv[n] = Wn.v - Wc.v;
                    dT[n] = temperature(Wn) - Tc;
                    ++n;
                };

                const int iL = (i == i_start_)     ? i_end_ - 1 : i - 1;
                const int iR = (i == i_end_ - 1)    ? i_start_   : i + 1;
                add(iL, j);
                add(iR, j);
                if (j > j_start_)     add(i, j - 1);
                if (j < j_end_ - 1)   add(i, j + 1);

                grad_u_[idx(i, j)] = least_squares_gradient(dx, dy, du, n);
                grad_v_[idx(i, j)] = least_squares_gradient(dx, dy, dv, n);
                grad_T_[idx(i, j)] = least_squares_gradient(dx, dy, dT, n);
            }
    }


    Conserved operator+(const Conserved& a, const Conserved& b) {
        return Conserved{
            a.rho  + b.rho,
            a.rhou + b.rhou,
            a.rhov + b.rhov,
            a.rhoE + b.rhoE
        };
    }

    Conserved operator-(const Conserved& a, const Conserved& b) {
        return Conserved{
            a.rho  - b.rho,
            a.rhou - b.rhou,
            a.rhov - b.rhov,
            a.rhoE - b.rhoE
        };
    }

    Conserved operator*(double s, const Conserved& a) {
        return Conserved{
            s * a.rho,
            s * a.rhou,
            s * a.rhov,
            s * a.rhoE
        };
    }

    Conserved operator*(const Conserved& a, double s) {
        return s * a;
    }

    // Residual Zeroing

    void RansSolver::zero_residual() { 

        for (auto& r : R_) { 

            r = Conserved{}; 

        }

    }

    void RansSolver::compute_interior_convective_residual() { 

        zero_residual(); 
        fill_ghost_cells(); 

        // i-direction faces: periodic around the o-grid

        for (int j = j_start_; j < j_end_; j++) { 

            for (int i = i_start_; i < i_end_; i++) { 

                const int iL = i; 

                const int iR_physical = (i == i_end_ - 1) ? i_start_ : i + 1; 

                const CellGeom& cL = cell_geom(iL, j); // left cell
                const CellGeom& cR = cell_geom(iR_physical, j); // right cell
                const FaceGeom& face = cL.face[1]; // i+ face of the left cell

                const Conserved& UL = state(iL, j); 
                const Conserved& UR = state(iR_physical, j); 

                const Conserved F = central_flux(UL, UR, face.nx, face.ny); 

                residual(iL, j) = residual(iL, j) + F * (face.length / cL.area); 
                residual(iR_physical, j) = residual(iR_physical, j) - F * (face.length / cR.area); 

            }
        }

        // j-direction interior faces: wall-normal direction, no periodicity

        for (int j = j_start_; j < j_end_ - 1; j++) { 

            for (int i = i_start_; i < i_end_; i++) { 

                const int jB = j; 
                const int jT = j + 1; 

                const CellGeom& cB = cell_geom(i, jB); // bottom cell
                const CellGeom& cT = cell_geom(i, jT); // top cell
                const FaceGeom& face = cB.face[2]; // j+ face of bottom cell

                const Conserved& UB = state(i, jB); 
                const Conserved& UT = state(i, jT); 

                const Conserved F = central_flux(UB, UT, face.nx, face.ny); 

                residual(i, jB) = residual(i, jB) + F * (face.length / cB.area); 
                residual(i, jT) = residual(i, jT) - F * (face.length / cT.area); 

            }
            
        }
 
    }

    double RansSolver::residual_linf_current() const { 

        double max_value = 0.0; 

        for (int j = j_start_; j< j_end_; j ++) { 

            for (int i = i_start_; i < i_end_; i++) { 

                const Conserved& r = residual(i, j); 

                max_value = std::max(max_value, std::fabs(r.rho));
                max_value = std::max(max_value, std::fabs(r.rhou));
                max_value = std::max(max_value, std::fabs(r.rhov));
                max_value = std::max(max_value, std::fabs(r.rhoE));

            }

        }
        
        return max_value; 

    }


    Primitive reflected_wall_state(const Primitive& Wi, double nx, double ny) { 

        if (!is_physical(Wi)) { 
            throw std::runtime_error("RANS: non-physical state in reflected wall state. "); 
        }

        const double un = Wi.u * nx + Wi.v * ny; 

        // ghost cells

        Primitive Wg; 
        Wg.rho = Wi.rho;
        Wg.p = Wi.p; 

        // speeds  
        Wg.u = Wi.u - 2.0 * un * nx; 
        Wg.v = Wi.v - 2.0 * un * ny; 

        return Wg; 

    }

    Primitive farfield_riemann_state(const Primitive& Wi, const Primitive& Winf, double nx, double ny) { 

        if (!is_physical(Wi)) { 
            throw std::runtime_error("RANS: Non-Physical interior state in farfield_riemann_state"); 
        }

        if (!is_physical(Winf)) { 
            throw std::runtime_error("RANS: Non-Physical freestream in farfield_riemann_state"); 
        }

        const double rhoi = Wi.rho; 
        const double ui = Wi.u; 
        const double vi = Wi.v; 
        const double pi = Wi.p; 

        const double rhoinf = Winf.rho; 
        const double uinf = Winf.u; 
        const double vinf = Winf.v; 
        const double pinf = Winf.p;

        const double ai = sound_speed(Wi); 
        const double ainf = sound_speed(Winf); 

        const double uni = ui * nx + vi * ny; 
        const double uninf = uinf * nx + vinf * ny;

        const double Vi = std::sqrt(ui * ui + vi * vi); 
        const double Mach_i = Vi / ai; 

        if (Mach_i >= 1.0) { 

            if (uni < 0.0) { 

                return Winf; 
                
            }

            return Wi; 
        }

        // subsonic inflow

        if (uni < 0.0) { 

            const double pb = 0.5 * (pinf + pi - rhoi * ai * (uninf - uni)); 
            const double rhob = rhoinf + (pb - pinf) / (ainf * ainf); 
            const double unb = uninf - (pinf - pb) / (rhoinf * ainf);
            const double ut_inf = -uinf * ny + vinf * nx; 

            Primitive Wb; 
            Wb.rho = rhob; 
            Wb.p = pb; 

            // speeds
            Wb.u = unb * nx - ut_inf * ny; 
            Wb.v = unb * ny + ut_inf * nx;

            return Wb; 

        }

        // subsonic outflow
        const double pb = pinf; 
        const double rhob = rhoi + (pb - pi) / (ai*ai); 
        const double unb = uni + (pi - pb) / (rhoi * ai); 
        const double ut_i = -ui * ny + vi * nx; 

        Primitive Wb; 
        Wb.p = pb; 
        Wb.rho = rhob; 

        // speeds 
        Wb.u = unb * nx - ut_i * ny; 
        Wb.v = unb * ny + ut_i * nx; 

        return Wb; 
    
    }

    void RansSolver::add_wall_convective_residual() { 

        const int j = j_start_; 

        for (int i = i_start_; i < i_end_; i++) { 

            const CellGeom& c = cell_geom(i, j); 
            const FaceGeom& face = c.face[0]; 

            const Conserved& Ui = state(i, j); 
            const Primitive Wi = conserved_to_primitive(Ui); 
            const Primitive Wg = reflected_wall_state(Wi, face.nx, face.ny); 

            Primitive Wface; 
            Wface.rho = 0.5 * (Wi.rho +  Wg.rho); 
            Wface.u = 0.5 * (Wi.u + Wg.u); 
            Wface.v = 0.5 * (Wi.v + Wg.v); 
            Wface.p = 0.5 * (Wi.p + Wg.p); 

            const Conserved Uface = primitive_to_conserved(Wface); 
            const Conserved F = normal_flux(Uface, face.nx, face.ny); 

            residual(i, j) = residual(i, j) + F * (face.length / c.area); 

        }

    }

    void RansSolver::add_farfield_convective_residual(const Primitive& Winf) { 

        const int j = j_end_ - 1; 

        for (int i = i_start_; i < i_end_; i++) { 

            const CellGeom& c = cell_geom(i, j); 
            const FaceGeom& face = c.face[2]; 

            const Conserved& Ui = state(i, j); 
            const Primitive Wi = conserved_to_primitive(Ui); 
            const Primitive Wb = farfield_riemann_state(Wi, Winf, face.nx, face.ny); 

            const Conserved F = central_flux(Ui, primitive_to_conserved(Wb), face.nx, face.ny); 

            residual(i, j) = residual(i, j) + F * (face.length / c.area); 

        }

    }

    void RansSolver::compute_full_convective_residual(const Primitive& Winf) { 

        compute_interior_convective_residual(); 
        add_wall_convective_residual(); 
        add_farfield_convective_residual(Winf); 

    }

    void RansSolver::add_interior_viscous_residual(double mu, double conductivity){ 

        if (mu < 0.0) { 

            throw std::runtime_error("RANS: Viscosity must be non-negative"); 
        }

        if (conductivity < 0.0) { 

            throw std::runtime_error("RANS: Conductivity must be non-negative"); 
        }

        compute_gradients(); 

        // i-direction faces: periodic around the O-grid. 

        for (int j =j_start_; j < j_end_; j++) { 

            for (int i = i_start_; i < i_end_; i++) { 

                const int iL = i; 
                const int iR = (i == i_end_ - 1 ) ? i_start_ : i + 1; 

                const CellGeom& cL = cell_geom(iL, j); 
                const CellGeom& cR = cell_geom(iR, j); 
                const FaceGeom& face = cL.face[1]; 

                const Primitive WL = conserved_to_primitive(state(iL, j)); 
                const Primitive WR = conserved_to_primitive(state(iR, j));

                // temperatures
                const double TL = temperature(WL); 
                const double TR = temperature(WR); 

                // gradients

                const double dx = cR.xc - cL.xc; 
                const double dy = cR.yc - cL.yc; 

                const Vec2 gu = corrected_face_gradient(grad_u(iL, j), grad_u(iR, j), WL.u, WR.u, dx, dy); 
                const Vec2 gv = corrected_face_gradient(grad_v(iL, j), grad_v(iR, j), WL.v, WR.v, dx, dy); 
                const Vec2 gT = corrected_face_gradient(grad_T(iL, j), grad_T(iR, j), TL, TR, dx, dy); 

                const double u_face = 0.5 * (WL.u + WR.u); 
                const double v_face = 0.5 * (WL.v + WR.v);

                const Conserved Fv = viscous_normal_flux(gu, gv, gT, u_face, v_face, mu, conductivity, face.nx, face.ny); 

                residual(iL, j) = residual(iL, j) - Fv * (face.length / cL.area); 
                residual(iR, j) = residual(iR, j) + Fv * (face.length / cR.area); 
                
            }
        }

        // j-direction interior faces: wall normal direction, no periodicitiy

        for (int j = j_start_; j < j_end_ - 1; j++) { 

            for (int i = i_start_; i < i_end_; i++) { 

                const int jB = j;
                const int jT = j + 1; 

                const CellGeom& cB = cell_geom(i, jB); 
                const CellGeom& cT = cell_geom(i, jT); 
                const FaceGeom& face = cB.face[2]; 

                const Primitive WB = conserved_to_primitive(state(i, jB)); 
                const Primitive WT = conserved_to_primitive(state(i, jT)); 

                // temperatures
                
                const double TB = temperature(WB); 
                const double TT = temperature(WT); 

                // gradients

                const double dx = cT.xc - cB.xc;
                const double dy = cT.yc - cB.yc; 

                const Vec2 gu = corrected_face_gradient(grad_u(i, jB), grad_u(i, jT), WB.u, WT.u, dx, dy); 
                const Vec2 gv = corrected_face_gradient(grad_v(i, jB), grad_v(i, jT), WB.v, WT.v, dx, dy);
                const Vec2 gT = corrected_face_gradient(grad_T(i, jB), grad_T(i, jT), TB, TT, dx, dy); 

                const double u_face = 0.5 * (WB.u + WT.u); 
                const double v_face = 0.5 * (WB.v + WT.v);

                const Conserved Fv = viscous_normal_flux(gu, gv, gT, u_face, v_face, mu, conductivity, face.nx, face.ny);
                residual(i, jB) = residual(i, jB) - Fv * (face.length / cB.area); 
                residual(i, jT) = residual(i, jT) + Fv * (face.length / cT.area); 

            }
        }
    }

    void RansSolver::add_wall_viscous_residual(double mu, double conductivity) { 

        if (mu < 0.0) { 

            throw std::runtime_error("RANS: Viscosity must be non-negative"); 

        }

        if (conductivity < 0.0) { 

            throw std::runtime_error("RANS: Conductivity must be non-negative"); 

        }

        const int j = j_start_; 

        for (int i = i_start_; i < i_end_; i++) { 

            const CellGeom& c = cell_geom(i, j); 
            const FaceGeom& face = c.face[0]; 

            const Primitive Wc = conserved_to_primitive(state(i, j)); 

            const double dx = c.xc - face.xc; 
            const double dy = c.yc - face.yc; 
            const double d = std::sqrt(dx*dx + dy*dy); 

            if (d <= 0.0) { 

                throw std::runtime_error("RANS: Zero wall distance in add_wall_viscous_residual"); 

            }

            /* 
            
            No slip wall: 
                u_wall = 0
                v_wall = 0

            Adiabatic wall: 

                dT/dn = 0
            
            Approximation: 
                Only the wall normal part of the velocity gradient is imposed here. Direction is froim wall face to cell center. 
            */

            const double ex = dx/d; 
            const double ey = dy/d; 

            const double du_dn = Wc.u / d; 
            const double dv_dn = Wc.v / d; 

            Vec2 grad_u_wall; 
            Vec2 grad_v_wall;
            Vec2 grad_T_wall; 

            grad_u_wall.x = du_dn * ex; 
            grad_u_wall.y = du_dn * ey; 

            grad_v_wall.x = dv_dn * ex; 
            grad_v_wall.y = dv_dn * ey; 

            grad_T_wall.x = 0.0; 
            grad_T_wall.y = 0.0; 

            const double u_face = 0.0; 
            const double v_face = 0.0;

            const Conserved Fv = viscous_normal_flux(grad_u_wall, grad_v_wall, grad_T_wall, u_face, v_face, mu, conductivity, face.nx, face.ny); 

            residual(i, j) = residual(i, j) - Fv * (face.length / c.area); 

        }

    }

    void RansSolver::compute_full_meanflow_residual(const Primitive& Winf, double mu, double conductivity) { 

        compute_full_convective_residual(Winf); 
        add_interior_viscous_residual(mu, conductivity); 
        add_wall_viscous_residual(mu, conductivity); 

    }

}
