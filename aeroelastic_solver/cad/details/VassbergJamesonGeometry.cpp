#include "VassbergJamesonGeometry.h"

#include <algorithm>
#include <cmath>

namespace
{
    constexpr double pi = 3.1415926535897932384626433832795;

    constexpr double t = 0.12;
    constexpr double x_te = 1.0089304115;
    constexpr double tau_te = 0.2818725;

    constexpr double z1_real = 1.0089304115;
    constexpr double z2_real = 0.0079337;

    constexpr double zeta1_real = 0.77043505;
    constexpr double zeta2_real = 0.24642903;
    constexpr double zeta_c_real = 0.4859156;

    double clamp_value(double value, double lo, double hi)
    {
        return std::max(lo, std::min(value, hi));
    }
}

namespace mesh
{
    double VassbergJamesonGeometry::reference_chord()
    {
        return 1.0;
    }

    double VassbergJamesonGeometry::trailing_edge_x()
    {
        return x_te;
    }

    double VassbergJamesonGeometry::thickness(double x)
    {
        x = clamp_value(x, 0.0, x_te);

        return 5.0 * t * (
            0.2969 * std::sqrt(x)
          - 0.1260 * x
          - 0.3516 * x * x
          + 0.2843 * x * x * x
          - 0.1015 * x * x * x * x
        );
    }

    VJPoint VassbergJamesonGeometry::surface_point(double s)
    {
        s = clamp_value(s, 0.0, 1.0);

        if (s <= 0.5)
        {
            const double q = s / 0.5;
            const double x = 0.5 * x_te * (1.0 + std::cos(pi * q));
            return VJPoint{x, thickness(x)};
        }

        const double q = (s - 0.5) / 0.5;
        const double x = 0.5 * x_te * (1.0 - std::cos(pi * q));
        return VJPoint{x, -thickness(x)};
    }

    std::complex<double> VassbergJamesonGeometry::kt_center()
    {
        return std::complex<double>(zeta_c_real, 0.0);
    }

    std::complex<double> VassbergJamesonGeometry::kt_forward(
        std::complex<double> z)
    {
        const std::complex<double> z1(z1_real, 0.0);
        const std::complex<double> z2(z2_real, 0.0);
        const std::complex<double> zeta1(zeta1_real, 0.0);
        const std::complex<double> zeta2(zeta2_real, 0.0);

        const double p = pi / (2.0 * pi - tau_te);

        const std::complex<double> a =
            std::pow((z - z1) / (z - z2), p);

        return (zeta1 - a * zeta2) / (1.0 - a);
    }

    std::complex<double> VassbergJamesonGeometry::kt_inverse(
        std::complex<double> zeta)
    {
        const std::complex<double> z1(z1_real, 0.0);
        const std::complex<double> z2(z2_real, 0.0);
        const std::complex<double> zeta1(zeta1_real, 0.0);
        const std::complex<double> zeta2(zeta2_real, 0.0);

        const double p = pi / (2.0 * pi - tau_te);

        const std::complex<double> w =
            std::pow((zeta - zeta1) / (zeta - zeta2), 1.0 / p);

        return (z1 - w * z2) / (1.0 - w);
    }
}