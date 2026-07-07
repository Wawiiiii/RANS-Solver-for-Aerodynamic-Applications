#include "euler_simulation.h"
#include "euler_multigrid_solver.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace aero {

    namespace {
        constexpr double PI = 3.14159265358979323846;

        // Strip leading/trailing whitespace.
        std::string trim(const std::string& s) {
            const char* ws = " \t\r\n";
            const size_t b = s.find_first_not_of(ws);
            if (b == std::string::npos) return "";
            return s.substr(b, s.find_last_not_of(ws) - b + 1);
        }

        // Typed parsers with key context for clear error messages.
        int as_int(const std::string& key, const std::string& v) {
            try { size_t p; const int x = std::stoi(v, &p); if (p == v.size()) return x; }
            catch (...) {}
            throw std::runtime_error("config: key '" + key + "' expects an integer, got '" + v + "'");
        }
        double as_double(const std::string& key, const std::string& v) {
            try { size_t p; const double x = std::stod(v, &p); if (p == v.size()) return x; }
            catch (...) {}
            throw std::runtime_error("config: key '" + key + "' expects a number, got '" + v + "'");
        }
        bool as_bool(const std::string& key, const std::string& v) {
            if (v == "true" || v == "1" || v == "yes" || v == "on")  return true;
            if (v == "false" || v == "0" || v == "no" || v == "off") return false;
            throw std::runtime_error("config: key '" + key + "' expects true/false, got '" + v + "'");
        }
    }

    SolverConfig load_solver_config(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) {
            throw std::runtime_error("Could not open config file: " + filename);
        }

        SolverConfig c;
        std::string line;
        int line_no = 0;

        while (std::getline(file, line)) {
            ++line_no;

            // Drop '#' comments, then trim.
            const size_t hash = line.find('#');
            if (hash != std::string::npos) line.erase(hash);
            line = trim(line);
            if (line.empty()) continue;

            const size_t eq = line.find('=');
            if (eq == std::string::npos) {
                throw std::runtime_error(
                    "config: line " + std::to_string(line_no) +
                    " is not 'key = value': " + line);
            }
            const std::string key = trim(line.substr(0, eq));
            const std::string val = trim(line.substr(eq + 1));

            if      (key == "mesh_file")        c.mesh_file = val;
            else if (key == "n_levels")         c.n_levels = as_int(key, val);
            else if (key == "mach")             c.mach = as_double(key, val);
            else if (key == "alpha_deg")        c.alpha_deg = as_double(key, val);
            else if (key == "cfl")              c.cfl = as_double(key, val);
            else if (key == "max_cycles")       c.max_cycles = as_int(key, val);
            else if (key == "print_every")      c.print_every = as_int(key, val);
            else if (key == "target_residual")  c.target_residual = as_double(key, val);
            else if (key == "k2")               c.k2 = as_double(key, val);
            else if (key == "k4")               c.k4 = as_double(key, val);
            else if (key == "n_pre")            c.n_pre = as_int(key, val);
            else if (key == "n_coarse")         c.n_coarse = as_int(key, val);
            else if (key == "n_post")           c.n_post = as_int(key, val);
            else if (key == "correction_omega") c.correction_omega = as_double(key, val);
            else if (key == "omega_decay")      c.omega_decay = as_double(key, val);
            else if (key == "coarse_cfl_decay") c.coarse_cfl_decay = as_double(key, val);
            else if (key == "bilinear_prolong") c.bilinear_prolong = as_bool(key, val);
            else if (key == "cycle_gamma")      c.cycle_gamma = as_int(key, val);
            else if (key == "use_fmg")          c.use_fmg = as_bool(key, val);
            else if (key == "n_fmg_cycles")     c.n_fmg_cycles = as_int(key, val);
            else if (key == "flowfield_file")   c.flowfield_file = val;
            else if (key == "surface_cp_file")  c.surface_cp_file = val;
            else throw std::runtime_error(
                "config: unknown key '" + key + "' at line " + std::to_string(line_no));
        }

        return c;
    }

    SolverResult EulerSimulation::solve() {
        const SolverConfig& c = config_;

        // Mesh + multigrid hierarchy.
        MeshGrid grid = load_indexed_mesh(c.mesh_file);
        EulerMultigridSolver mg(grid.x, grid.y, c.n_levels);

        // Freestream initial condition.
        const double alpha_rad = c.alpha_deg * PI / 180.0;
        const Primitive Winf = freestream_state(c.mach, alpha_rad);

        std::cout << std::scientific << std::setprecision(12)
                  << "\nFreestream (M=" << c.mach << ", alpha=" << c.alpha_deg << " deg):\n"
                  << "  rho=" << Winf.rho << "  u=" << Winf.u
                  << "  v=" << Winf.v << "  p=" << Winf.p
                  << "  a=" << sound_speed(Winf) << "\n";

        mg.set_uniform_state(Winf);

        JSTParameters jst;
        jst.k2 = c.k2;
        jst.k4 = c.k4;

        // Drive to convergence (FMG start-up or plain cycles).
        const int cycles = c.use_fmg
            ? mg.run_full_multigrid(
                  Winf, jst, c.cfl, c.max_cycles, c.print_every, c.target_residual,
                  c.n_pre, c.n_coarse, c.n_post, c.n_fmg_cycles,
                  c.correction_omega, c.omega_decay, c.coarse_cfl_decay,
                  c.bilinear_prolong, c.cycle_gamma)
            : mg.run_v_cycles(
                  Winf, jst, c.cfl, c.max_cycles, c.print_every, c.target_residual,
                  c.n_pre, c.n_coarse, c.n_post,
                  c.correction_omega, c.omega_decay, c.coarse_cfl_decay,
                  c.bilinear_prolong, c.cycle_gamma);

        // Outputs.
        if (!c.flowfield_file.empty())  mg.fine().write_flowfield(c.flowfield_file);
        if (!c.surface_cp_file.empty()) mg.fine().write_surface_cp(c.surface_cp_file, Winf);

        // Result.
        mg.fine().compute_full_residual_jst(Winf, jst);
        SolverResult res;
        res.coeffs         = mg.fine().compute_aero_coefficients(Winf, alpha_rad);
        res.cycles         = cycles;
        res.final_residual = mg.fine().residual_linf_current();
        res.converged      = (cycles < c.max_cycles);
        return res;
    }
}
