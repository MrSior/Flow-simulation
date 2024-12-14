
#ifndef FLUID_PROJECT_VECTORFIELD_H
#define FLUID_PROJECT_VECTORFIELD_H

#include <utility>
#include <ranges>
#include <cassert>

#include "common.h"
#include "customArray.h"

namespace Physics {
    template<typename T, int N, int M>
    struct VectorField {
        Array<std::array<T, deltas.size()>, N, M> v;

        T &add(int x, int y, int dx, int dy, T dv) {
            return get(x, y, dx, dy) += dv;
        }

        T &get(int x, int y, int dx, int dy) {
            size_t i = std::ranges::find(deltas, std::pair(dx, dy)) - deltas.begin();
            assert(i < deltas.size());
            return v[x][y][i];
        }
    };
}

#endif //FLUID_PROJECT_VECTORFIELD_H
