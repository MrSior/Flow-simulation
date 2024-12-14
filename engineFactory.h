
#ifndef FLUID_PROJECT_ENGINEFACTORY_H
#define FLUID_PROJECT_ENGINEFACTORY_H

#include <cinttypes>
#include <type_traits>
#include <array>
#include <iostream>
#include "FluidEngine.h"

//#define TYPES FLOAT,DOUBLE,FIXED(32, 9),FAST_FIXED(40, 5),FAST_FIXED(50, 16)

#define FLOAT 1
#define DOUBLE 2
#define FIXED(n, k) ((n) * 1000 + (k))
#define FAST_FIXED(n, k) ((n) * 100000 + (k))
#define PAIR(n, k) std::pair<int, int>((n), (k))
#define S(n, m) std::pair<int, int>((n), (m))
//#define SIZES S(35, 286),S(36, 84)


#ifndef TYPES
#error "Not defined any types"
#endif

#ifndef SIZES
#define SIZES S(-1, -1)
#endif

namespace Physics {
    template<int n>
    struct TypeHolder {
        using type = std::conditional_t<n == 1, float,
                std::conditional_t<n == 2, double,
                        std::conditional_t<(n > 100000), Fixed<n / 100000, n % 100000, true>,
                                std::conditional_t<(n > 1000), Fixed<n / 1000, n % 1000>, void>>>>;
    };

    template<int n>
    using getType = TypeHolder<n>::type;


    constexpr auto GenerateAllCombo() {
        constexpr std::pair<int, int> givenSizes[] = {SIZES, {-1, -1}};
        constexpr int sizesCnt = sizeof(givenSizes) / sizeof(std::pair<int, int>);

        constexpr auto givenTypes = std::array{TYPES};
        constexpr std::array<std::pair<int, int>, sizesCnt> sizes = {std::pair<int, int>(-1, -1), SIZES};
        std::array<std::tuple<int, int, int, int, int>,
                givenTypes.size() * givenTypes.size() * givenTypes.size() * sizesCnt> res = {};

        int idx = 0;
        for (int pType: givenTypes) {
            for (int vType: givenTypes) {
                for (int vfType: givenTypes) {
                    for (auto field: sizes) {
                        res[idx] = {pType, vType, vfType, field.first, field.second};
                        ++idx;
                    }
                }
            }
        }
        return res;
    }

    constexpr auto allCombos = GenerateAllCombo();
    std::array<std::shared_ptr<Physics::IEngine>(*)(), allCombos.size()> generateEngine;


    template<int idx>
    struct EngineGeneratorIdx {
        EngineGeneratorIdx<idx - 1> x;

        EngineGeneratorIdx() {
            generateEngine[idx - 1] = generate;
        }

        static std::shared_ptr<Physics::IEngine> generate() {
            return std::make_shared<Physics::Engine<getType<std::get<0>(allCombos[idx - 1])>,
                    getType<std::get<1>(allCombos[idx - 1])>,
                    getType<std::get<2>(allCombos[idx - 1])>,
                    std::get<3>(allCombos[idx - 1]),
                    std::get<4>(allCombos[idx - 1])>>();
        }
    };

    template<>
    struct EngineGeneratorIdx<0> {
        EngineGeneratorIdx() = default;
    };

    [[maybe_unused]] EngineGeneratorIdx<allCombos.size()> generator{};
}

std::shared_ptr<Physics::IEngine> ProduceEngine(int pType, int vType, int vfType, int n, int m) {
    auto itr = std::find(Physics::allCombos.begin(),
                         Physics::allCombos.end(),
                         std::tuple(pType, vType, vfType, n, m));
    std::shared_ptr<Physics::IEngine> engine;
    if (itr != Physics::allCombos.end()) {
        engine = Physics::generateEngine[itr - Physics::allCombos.begin()]();
    } else {
        itr = std::find(Physics::allCombos.begin(),
                             Physics::allCombos.end(),
                             std::tuple(pType, vType, vfType, -1, -1));
        if (itr == Physics::allCombos.end()) {
            throw std::invalid_argument("unknown types");
        }
        engine = Physics::generateEngine[itr - Physics::allCombos.begin()]();
    }
    return engine;
}

#endif //FLUID_PROJECT_ENGINEFACTORY_H
