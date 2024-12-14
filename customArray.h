//
// Created by Семён on 11.12.2024.
//

#ifndef FLUID_PROJECT_CUSTOMARRAY_H
#define FLUID_PROJECT_CUSTOMARRAY_H

#include <vector>
#include <cstring>

namespace Physics {
    template<typename T, int valN, int valM>
    struct Array{
        T arr[valN][valM]{};
        int N = valN;
        int M = valM;

        void init(int n, int m) {}

        T *operator[](int n) {
            return arr[n];
        }

        Array &operator=(const Array &other) {
            if (this == &other) {
                return *this;
            }
            std::memcpy(arr, other.arr, sizeof(arr));
            return *this;
        }
    };

    template<typename T>
    struct Array<T, -1, -1>{
        std::vector<std::vector<T>> arr{};
        int N = 0;
        int M = 0;

        void init(int n, int m) {
            N = n;
            M = m;
            arr.assign(n, std::vector<T>(m, T{}));
        }

        std::vector<T> &operator[](int n) {
            return arr[n];
        }

        Array &operator=(const Array &other) {
            if (this == &other) {
                return *this;
            }
            arr = other.arr;
            return *this;
        }
    };

}


#endif //FLUID_PROJECT_CUSTOMARRAY_H
