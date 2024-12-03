#include "bits/stdc++.h"
using namespace std;

//constexpr size_t N = 36, M = 84;
 constexpr size_t N = 14, M = 5;
constexpr size_t T = 1'000'000;
constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

 char field[N][M + 1] = {
     "#####",
     "#.  #",
     "#.# #",
     "#.# #",
     "#.# #",
     "#.# #",
     "#.# #",
     "#.# #",
     "#...#",
     "#####",
     "#   #",
     "#   #",
     "#   #",
     "#####",
 };

//char field[N][M + 1] = {
//        "####################################################################################",
//        "#                                                                                  #",
//        "#                                                                                  #",
//        "#                                                                                  #",
//        "#                                                                                  #",
//        "#                                                                                  #",
//        "#                                       .........                                  #",
//        "#..............#            #           .........                                  #",
//        "#..............#            #           .........                                  #",
//        "#..............#            #           .........                                  #",
//        "#..............#            #                                                      #",
//        "#..............#            #                                                      #",
//        "#..............#            #                                                      #",
//        "#..............#            #                                                      #",
//        "#..............#............#                                                      #",
//        "#..............#............#                                                      #",
//        "#..............#............#                                                      #",
//        "#..............#............#                                                      #",
//        "#..............#............#                                                      #",
//        "#..............#............#                                                      #",
//        "#..............#............#                                                      #",
//        "#..............#............#                                                      #",
//        "#..............#............################                     #                 #",
//        "#...........................#....................................#                 #",
//        "#...........................#....................................#                 #",
//        "#...........................#....................................#                 #",
//        "##################################################################                 #",
//        "#                                                                                  #",
//        "#                                                                                  #",
//        "#                                                                                  #",
//        "#                                                                                  #",
//        "#                                                                                  #",
//        "#                                                                                  #",
//        "#                                                                                  #",
//        "#                                                                                  #",
//        "####################################################################################",
//};

template<uint8_t N, uint8_t K>
struct Fixed {
    constexpr static uint8_t n = N;
    constexpr static uint8_t k = K;

    using val_t =   std::conditional_t<n <= 8, int8_t,
                    std::conditional_t<n <= 16, int16_t,
                    std::conditional_t<n <= 32, int32_t,
                    std::conditional_t<n <= 64, int64_t, void>>>>;

    static_assert(!std::is_same_v<val_t, void>, "n must be less than 64");

    constexpr Fixed(val_t v): v(v << k) {}
    constexpr Fixed(float f): v(f * (1 << k)) {}
    constexpr Fixed(double f): v(f * (1 << k)) {}
    constexpr Fixed(): v(0) {}

    static constexpr Fixed from_raw(val_t x) {
        Fixed ret;
        ret.v = x;
        return ret;
    }

    val_t v;

    auto operator<=>(const Fixed&) const = default;
    bool operator==(const Fixed&) const = default;

    template<uint8_t otherN, uint8_t otherK>
    auto operator+(Fixed<otherN, otherK> other) {
        if constexpr (n > other.n) {
            return Fixed::from_raw(this->v + other.v);
        } else {
            return decltype(other)::from_raw(this->v + other.v);
        }
    }

    template<uint8_t otherN, uint8_t otherK>
    auto operator-(Fixed<otherN, otherK> other) {
        if constexpr (n > other.n) {
            return decltype(*this)::from_raw(v - other.v);
        } else {
            return decltype(other)::from_raw(v - other.v);
        }
    }

    template<uint8_t otherN, uint8_t otherK>
    auto operator*(Fixed<otherN, otherK> other) {
        if constexpr (n > other.n) {
            return decltype(*this)::from_raw(((int64_t) v * other.v) >> 16);
        } else {
            return decltype(other)::from_raw(((int64_t) v * other.v) >> 16);
        }
    }

    template<uint8_t otherN, uint8_t otherK>
    auto operator/(Fixed<otherN, otherK> other) {
        if constexpr (n > other.n) {
            return decltype(*this)::from_raw(((int64_t) v << 16) / other.v);
        } else {
            return decltype(other)::from_raw(((int64_t) v << 16) / other.v);
        }
    }

    auto operator/(int other_v) {
        if constexpr (n > 32) {
            return decltype(*this)::from_raw(((int64_t) v << 16) / other_v);
        } else {
            return Fixed<32, 16>::from_raw(((int64_t) v << 16) / other_v);
        }
    }

    template<uint8_t otherN, uint8_t otherK>
    Fixed<n, k>& operator+=(Fixed<otherN, otherK> other) {
        return *this = *this + other;
    }

    template<uint8_t otherN, uint8_t otherK>
    Fixed<n, k>& operator-=(Fixed<otherN, otherK> other) {
        return *this = *this - other;
    }

    template<uint8_t otherN, uint8_t otherK>
    Fixed<n, k>& operator*=(Fixed<otherN, otherK> other) {
        return *this = *this * other;
    }

    Fixed<n, k>& operator*=(float other) {
        return *this = *this * Fixed<32,16>(other);
    }

    template<uint8_t otherN, uint8_t otherK>
    Fixed<n, k>& operator/=(Fixed<otherN, otherK> other) {
        return *this = *this / other;
    }

    Fixed<n, k> operator-(Fixed<n, k> other) {
        return Fixed<n, k>::from_raw(-v);
    }

    ostream &operator<<(ostream &out) {
        return out << v / (double) (1 << 16);
    }
};

//Fixed abs(Fixed x) {
//    if (x.v < 0) {
//        x.v = -x.v;
//    }
//    return x;
//}

Fixed<32, 16> rho[256];

Fixed<32, 16> p[N][M]{}, old_p[N][M];

struct VectorField {
    array<Fixed<32, 16>, deltas.size()> v[N][M];
    Fixed<32, 16> &add(int x, int y, int dx, int dy, Fixed<32, 16> dv) {
        return get(x, y, dx, dy) += dv;
    }

    Fixed<32, 16> &get(int x, int y, int dx, int dy) {
        size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
        assert(i < deltas.size());
        return v[x][y][i];
    }
};

VectorField velocity{}, velocity_flow{};
int last_use[N][M]{};
int UT = 0;


mt19937 rnd(1337);

tuple<Fixed<32, 16>, bool, pair<int, int>> propagate_flow(int x, int y, Fixed<32, 16> lim) {
    last_use[x][y] = UT - 1;
    Fixed<32, 16> ret = 0;
    for (auto [dx, dy] : deltas) {
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
            auto cap = velocity.get(x, y, dx, dy);
            auto flow = velocity_flow.get(x, y, dx, dy);
            if (flow == cap) {
                continue;
            }
            // assert(v >= velocity_flow.get(x, y, dx, dy));
            auto vp = min(lim, cap - flow);
            if (last_use[nx][ny] == UT - 1) {
                velocity_flow.add(x, y, dx, dy, vp);
                last_use[x][y] = UT;
                // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
                return {vp, 1, {nx, ny}};
            }
            auto [t, prop, end] = propagate_flow(nx, ny, vp);
            ret += t;
            if (prop) {
                velocity_flow.add(x, y, dx, dy, t);
                last_use[x][y] = UT;
                // cerr << x << " " << y << " -> " << nx << " " << ny << " " << t << " / " << lim << "\n";
                return {t, prop && end != pair(x, y), end};
            }
        }
    }
    last_use[x][y] = UT;
    return {ret, 0, {0, 0}};
}

Fixed<32, 16> random01() {
    return Fixed<32, 16>::from_raw((rnd() & ((1 << 16) - 1)));
}

void propagate_stop(int x, int y, bool force = false) {
    if (!force) {
        bool stop = true;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > 0) {
                stop = false;
                break;
            }
        }
        if (!stop) {
            return;
        }
    }
    last_use[x][y] = UT;
    for (auto [dx, dy] : deltas) {
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > 0) {
            continue;
        }
        propagate_stop(nx, ny);
    }
}

Fixed<32, 16> move_prob(int x, int y) {
    Fixed<32, 16> sum = 0;
    for (size_t i = 0; i < deltas.size(); ++i) {
        auto [dx, dy] = deltas[i];
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
            continue;
        }
        auto v = velocity.get(x, y, dx, dy);
        if (v < 0) {
            continue;
        }
        sum += v;
    }
    return sum;
}

struct ParticleParams {
    char type;
    Fixed<32, 16> cur_p;
    array<Fixed<32, 16>, deltas.size()> v;

    void swap_with(int x, int y) {
        swap(field[x][y], type);
        swap(p[x][y], cur_p);
        swap(velocity.v[x][y], v);
    }
};

bool propagate_move(int x, int y, bool is_first) {
    last_use[x][y] = UT - is_first;
    bool ret = false;
    int nx = -1, ny = -1;
    do {
        std::array<Fixed<32, 16>, deltas.size()> tres;
        Fixed<32, 16> sum = 0;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                tres[i] = sum;
                continue;
            }
            auto v = velocity.get(x, y, dx, dy);
            if (v < 0) {
                tres[i] = sum;
                continue;
            }
            sum += v;
            tres[i] = sum;
        }

        if (sum == 0) {
            break;
        }

        Fixed<32, 16> p = random01() * sum;
        size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

        auto [dx, dy] = deltas[d];
        nx = x + dx;
        ny = y + dy;
        assert(velocity.get(x, y, dx, dy) > 0 && field[nx][ny] != '#' && last_use[nx][ny] < UT);

        ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
    } while (!ret);
    last_use[x][y] = UT;
    for (size_t i = 0; i < deltas.size(); ++i) {
        auto [dx, dy] = deltas[i];
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < 0) {
            propagate_stop(nx, ny);
        }
    }
    if (ret) {
        if (!is_first) {
            ParticleParams pp{};
            pp.swap_with(x, y);
            pp.swap_with(nx, ny);
            pp.swap_with(x, y);
        }
    }
    return ret;
}

int dirs[N][M]{};

int main() {
    rho[' '] = 0.01;
    rho['.'] = 1000;
    Fixed<32, 16> g = 0.1;

    for (size_t x = 0; x < N; ++x) {
        for (size_t y = 0; y < M; ++y) {
            if (field[x][y] == '#')
                continue;
            for (auto [dx, dy] : deltas) {
                dirs[x][y] += (field[x + dx][y + dy] != '#');
            }
        }
    }

    for (size_t i = 0; i < T; ++i) {

        Fixed<32, 16> total_delta_p = 0;
        // Apply external forces
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                if (field[x + 1][y] != '#')
                    velocity.add(x, y, 1, 0, g);
            }
        }

        // Apply forces from p
        memcpy(old_p, p, sizeof(p));
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    int nx = x + dx, ny = y + dy;
                    if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
                        auto delta_p = old_p[x][y] - old_p[nx][ny];
                        auto force = delta_p;
                        auto &contr = velocity.get(nx, ny, -dx, -dy);
                        if (contr * rho[(int) field[nx][ny]] >= force) {
                            contr -= force / rho[(int) field[nx][ny]];
                            continue;
                        }
                        force -= contr * rho[(int) field[nx][ny]];
                        contr = 0;
                        velocity.add(x, y, dx, dy, force / rho[(int) field[x][y]]);
                        p[x][y] -= force / dirs[x][y];
                        total_delta_p -= force / dirs[x][y];
                    }
                }
            }
        }

        // Make flow from velocities
        velocity_flow = {};
        bool prop = false;
        do {
            UT += 2;
            prop = 0;
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] != '#' && last_use[x][y] != UT) {
                        auto [t, local_prop, _] = propagate_flow(x, y, 1);
                        if (t > 0) {
                            prop = 1;
                        }
                    }
                }
            }
        } while (prop);

        // Recalculate p with kinetic energy
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    auto old_v = velocity.get(x, y, dx, dy);
                    auto new_v = velocity_flow.get(x, y, dx, dy);
                    if (old_v > 0) {
                        assert(new_v <= old_v);
                        velocity.get(x, y, dx, dy) = new_v;
                        auto force = (old_v - new_v) * rho[(int) field[x][y]];
                        if (field[x][y] == '.')
                            force *= 0.8;
                        if (field[x + dx][y + dy] == '#') {
                            p[x][y] += force / dirs[x][y];
                            total_delta_p += force / dirs[x][y];
                        } else {
                            p[x + dx][y + dy] += force / dirs[x + dx][y + dy];
                            total_delta_p += force / dirs[x + dx][y + dy];
                        }
                    }
                }
            }
        }

        UT += 2;
        prop = false;
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] != '#' && last_use[x][y] != UT) {
                    if (random01() < move_prob(x, y)) {
                        prop = true;
                        propagate_move(x, y, true);
                    } else {
                        propagate_stop(x, y, true);
                    }
                }
            }
        }

        prop = true;
        if (prop) {
            cout << "Tick " << i << ":\n";
            for (size_t x = 0; x < N; ++x) {
                cout << field[x] << "\n";
            }
        }

    }
}