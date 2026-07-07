#pragma once

#include <complex>

namespace mesh
{
    struct VJPoint
    {
        double x = 0.0;
        double y = 0.0;
    };

    class VassbergJamesonGeometry
    {
    public:
        static double reference_chord();
        static double trailing_edge_x();

        static double thickness(double x);
        static VJPoint surface_point(double s);

        static std::complex<double> kt_center();
        static std::complex<double> kt_forward(std::complex<double> z);
        static std::complex<double> kt_inverse(std::complex<double> zeta);
    };
}