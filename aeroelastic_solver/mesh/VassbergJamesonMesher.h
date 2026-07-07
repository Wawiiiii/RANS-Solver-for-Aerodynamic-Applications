#pragma once

#include "euler_mesher.h"

#include <cstddef>

namespace mesh
{
    class VassbergJamesonMesher
    {
    public:
        static Mesh2D generate(size_t n_cells);
    };
}