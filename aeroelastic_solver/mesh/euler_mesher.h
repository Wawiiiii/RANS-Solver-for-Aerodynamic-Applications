#pragma once

#include <vector>
#include <cmath>
#include <array>

namespace mesh
{
    struct Mesh2D
    {
        std::vector<double> x_grid;
        std::vector<double> y_grid;
        size_t n_t, n_r;
    };

    struct Profile
    {
        std::vector<double> x_points;
        std::vector<double> y_points;
    };

    class EulerMesher
    {
    public:
        // initialize the solver from a mesh or a profile
        EulerMesher(Mesh2D& mesh);
        EulerMesher(const Profile& profile, size_t n_t, size_t n_r, double farfield);
        Mesh2D get_mesh() const;

        // can build a mesh from a profile
        static Mesh2D initialize_mesh(const Profile& profile, size_t n_t, size_t n_r, double farfield);

        // apply eliptical smoothing to the mesh until convergence
        void apply_elliptical_smoothing();

    private:
        static constexpr double tol = 1e-10;
        static constexpr size_t max_iter = 1e5;

        std::vector<double> x_grid;
        std::vector<double> y_grid;

        std::vector<double> x_intergrid;
        std::vector<double> y_intergrid;

        size_t n_t, n_r;
        inline size_t idx(int i, int j) {
            return (size_t) ( i + j * n_t );
        };

        double elliptical_smoothing_step();
    };
};
