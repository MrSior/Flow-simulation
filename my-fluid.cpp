#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <tuple>
#include <algorithm>
#include <sstream>

#include "Fixed.h"

using namespace std;

//constexpr size_t N = 14, M = 5;
constexpr size_t N = 36, M = 84;

constexpr size_t T = 1'000'000;
constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

//char field[N][M + 1] = {
//        "#####",
//        "#.  #",
//        "#.# #",
//        "#.# #",
//        "#.# #",
//        "#.# #",
//        "#.# #",
//        "#.# #",
//        "#...#",
//        "#####",
//        "#   #",
//        "#   #",
//        "#   #",
//        "#####",
//};

char field[N][M + 1] = {
    "####################################################################################",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                       .........                                  #",
    "#..............#            #           .........                                  #",
    "#..............#            #           .........                                  #",
    "#..............#            #           .........                                  #",
    "#..............#            #                                                      #",
    "#..............#            #                                                      #",
    "#..............#            #                                                      #",
    "#..............#            #                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............################                     #                 #",
    "#...........................#....................................#                 #",
    "#...........................#....................................#                 #",
    "#...........................#....................................#                 #",
    "##################################################################                 #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "####################################################################################",
};

//template <uint8_t N, uint8_t K>
//class Fixed {
//    static_assert(N > K, "N must be greater than K");
//    using value_t = std::conditional_t<
//            N <= 8, int8_t,
//            std::conditional_t<
//                    N <= 16, int16_t,
//                    std::conditional_t<
//                            N <= 32, int32_t,
//                            std::conditional_t<N == 64, int64_t, void>>>>;
//
//
//    static constexpr int32_t scale = 1 << K;
//
//    value_t v;
//
//public:
//    constexpr Fixed(int v = 0) : v(v << K) {}
//    constexpr Fixed(float f) : v(static_cast<value_t>(f * scale)) {}
//    constexpr Fixed(double f) : v(static_cast<value_t>(f * scale)) {}
//
//    static constexpr Fixed from_raw(value_t x) {
//        Fixed ret;
//        ret.v = x;
//        return ret;
//    }
//
//    constexpr auto operator<=>(const Fixed&) const = default;
//    constexpr bool operator==(const Fixed&) const = default;
//
//    Fixed operator+(Fixed other) const { return from_raw(v + other.v); }
//    Fixed operator-(Fixed other) const { return from_raw(v - other.v); }
//    Fixed operator*(Fixed other) const {
//        return from_raw(static_cast<int64_t>(v) * other.v / scale);
//    }
//    Fixed operator/(Fixed other) const {
//        return from_raw(static_cast<int64_t>(v) * scale / other.v);
//    }
//
//    Fixed& operator+=(Fixed other) { return *this = *this + other; }
//    Fixed& operator-=(Fixed other) { return *this = *this - other; }
//    Fixed& operator*=(Fixed other) { return *this = *this * other; }
//    Fixed& operator/=(Fixed other) { return *this = *this / other; }
//
//    Fixed operator-() const { return from_raw(-v); }
//
//    constexpr Fixed abs() const { return from_raw(v < 0 ? -v : v); }
//
//    friend ostream& operator<<(ostream& out, const Fixed& x) {
//        return out << static_cast<double>(x.v) / scale;
//    }
//};






template <typename pType, typename velocityType, typename velocityFlowtype>
class EngineImpl{
public:
    using p_t = pType;
    using vel_t = velocityType;
    using vflow_t = velocityFlowtype;

//    p_t rho[256];
    p_t p[N][M]{}, old_p[N][M];

    template<typename T>
    struct VectorField {
        array<T, deltas.size()> v[N][M];
        T& add(int x, int y, int dx, int dy, T dv) {
            return get(x, y, dx, dy) += dv;
        }

        T& get(int x, int y, int dx, int dy) {
            size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
            assert(i < deltas.size());
            return v[x][y][i];
        }
    };

    VectorField<vel_t> velocity{};
    VectorField<vflow_t> velocity_flow{};
    int last_use[N][M]{};
    int UT = 0;

    mt19937 rnd = mt19937(1337);

    tuple<vflow_t, bool, pair<int, int>> propagate_flow(int x, int y, vflow_t lim) {
        last_use[x][y] = UT - 1;
        vflow_t ret = 0;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
                auto cap = velocity.get(x, y, dx, dy);
                auto flow = velocity_flow.get(x, y, dx, dy);
                if (flow == vflow_t(cap)) {
                    continue;
                }
                // assert(v >= velocity_flow.get(x, y, dx, dy));
                auto vp = min(lim, vflow_t(cap - flow));
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

    p_t random01() {
        if constexpr (std::is_same_v<p_t, float> || std::is_same_v<p_t, double>) {
            return rnd() & ((1 << 16) - 1);
        } else {
            return p_t::from_raw((rnd() & ((1 << p_t::k) - 1)));
        }
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

    p_t move_prob(int x, int y) {
        p_t sum = 0;
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
        p_t cur_p;
        array<vel_t , deltas.size()> v;

        void swap_with(int x, int y, VectorField<vel_t>& vel, p_t p_ptr[N][M]) {
            swap(field[x][y], type);
            swap(p_ptr[x][y], cur_p);
            swap(vel.v[x][y], v);
        }
    };

    bool propagate_move(int x, int y, bool is_first) {
        last_use[x][y] = UT - is_first;
        bool ret = false;
        int nx = -1, ny = -1;
        do {
            std::array<p_t, deltas.size()> tres;
            p_t sum = 0;
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

            p_t p = random01() * sum;
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
                pp.swap_with(x, y, velocity, p);
                pp.swap_with(nx, ny, velocity, p);
                pp.swap_with(x, y, velocity, p);
            }
        }
        return ret;
    }

    int dirs[N][M]{};



    void run(Fixed<32, 16> rho[256], Fixed<32, 16>& g) {

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

            p_t total_delta_p = 0;
            // Apply external forces
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    if (field[x + 1][y] != '#')
                        velocity.add(x, y, 1, 0, vel_t(g));
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
                            if (contr * rho[(int) field[nx][ny]] >= vel_t(force)) {
                                contr -= force / rho[(int) field[nx][ny]];
                                continue;
                            }
                            force -= contr * rho[(int) field[nx][ny]];
                            contr = 0;
                            velocity.add(x, y, dx, dy, vel_t(force / rho[(int) field[x][y]]));
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
                            assert(new_v <= vflow_t(old_v));
                            velocity.get(x, y, dx, dy) = vel_t(new_v);
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

            if (prop) {
                cout << "Tick " << i << ":\n";
                for (size_t x = 0; x < N; ++x) {
                    cout << field[x] << "\n";
                }
            }
        }
    }

    int curTick = 0;
    Fixed<32, 16> rho[256];
    Fixed<32, 16> g;
    void initField(Fixed<32, 16> rho[256], Fixed<32, 16>& g) {
        memcpy(this->rho, rho, 256 * sizeof (rho[0]));
        this->g = g;
        curTick = 0;
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    dirs[x][y] += (field[x + dx][y + dy] != '#');
                }
            }
        }
    }

    std::string next() {
        p_t total_delta_p = 0;
        // Apply external forces
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                if (field[x + 1][y] != '#')
                    velocity.add(x, y, 1, 0, vel_t(g));
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
                        if (contr * rho[(int) field[nx][ny]] >= vel_t(force)) {
                            contr -= force / rho[(int) field[nx][ny]];
                            continue;
                        }
                        force -= contr * rho[(int) field[nx][ny]];
                        contr = 0;
                        velocity.add(x, y, dx, dy, vel_t(force / rho[(int) field[x][y]]));
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
                        assert(new_v <= vflow_t(old_v));
                        velocity.get(x, y, dx, dy) = vel_t(new_v);
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

        std::stringstream ss("");
        if (prop) {
            ss << "Tick " << curTick++ << ":\n";
            for (size_t x = 0; x < N; ++x) {
                ss << field[x] << "\n";
            }
        }
        return ss.str();
    }
};




int main() {
    Fixed<32, 16> rho[256];
    rho[' '] = Fixed<32, 16>(0.01);
    rho['.'] = Fixed<32, 16>(1000);
    auto g = Fixed<32, 16>(0.1);

    EngineImpl<Fixed<32, 9>, Fixed<32, 9>, double> engine;
//    engine.run(rho, g);
    engine.initField(rho, g);
    for (int i = 0; i < 1000000; ++i) {
        std::cout << engine.next();
    }
}