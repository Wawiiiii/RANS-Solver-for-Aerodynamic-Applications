#include "RansMesher.h"
#include "rans_solver.h"

#include <cmath>
#include <iostream>

int main()
{
    mesh::RansMeshParams params;
    params.n_cells = 128;
    params.n_radial = 96;
    params.reynolds = 1.0e6;
    params.y_plus = 1.0;
    params.growth = 1.18;
    params.farfield = 20.0;

    const mesh::Mesh2D mesh = mesh::RansMesher::generate(params);

    rans::RansSolver solver(
        mesh.x_grid,
        mesh.y_grid,
        static_cast<int>(mesh.n_t),
        static_cast<int>(mesh.n_r));

    const double mach = 0.20;
    const double rho_inf = 1.0;
    const double p_inf = 1.0 / rans::GAMMA;
    const double a_inf = std::sqrt(rans::GAMMA * p_inf / rho_inf);
    const double u_inf = mach * a_inf;
    const double v_inf = 0.0;

    const rans::Primitive Winf{rho_inf, u_inf, v_inf, p_inf};

    const double chord = 1.0;
    const double mu = rho_inf * u_inf * chord / params.reynolds;
    const double prandtl = 0.72;
    const double cp = rans::GAMMA / (rans::GAMMA - 1.0);
    const double conductivity = mu * cp / prandtl;

    const double cfl = 1.0e-5;
    const int max_iterations = 20;
    const int print_every = 1;
    const double target_residual = 0.0;

    std::cout << "RANS run setup\n";
    std::cout << "  cells       : " << solver.ni_cells() << " x " << solver.nj_cells() << "\n";
    std::cout << "  Mach        : " << mach << "\n";
    std::cout << "  Reynolds    : " << params.reynolds << "\n";
    std::cout << "  mu          : " << mu << "\n";
    std::cout << "  conductivity: " << conductivity << "\n";
    std::cout << "  CFL         : " << cfl << "\n";

    solver.set_uniform_state(Winf);

    const int final_iter = solver.run_pseudo_time_iterations(
        Winf,
        mu,
        conductivity,
        cfl,
        max_iterations,
        print_every,
        target_residual);

    solver.write_flowfield("rans_flowfield.dat");
    solver.write_wall_data("rans_wall.dat", Winf, mu);

    std::cout << "Finished after " << final_iter << " pseudo-time iterations\n";
    std::cout << "Final residual_linf = " << solver.residual_linf_current() << "\n";
    std::cout << "Wrote rans_flowfield.dat and rans_wall.dat\n";

    return 0;
}
