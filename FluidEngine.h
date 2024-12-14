
#ifndef FLUID_PROJECT_FLUIDENGINE_H
#define FLUID_PROJECT_FLUIDENGINE_H

#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <tuple>
#include <algorithm>
#include <sstream>
#include <optional>
#include <fstream>

#include "vectorField.h"
#include "customArray.h"
#include "common.h"
#include "Fixed.h"

using namespace std;


namespace Physics {
    class IEngine {
    public:
        virtual void next(std::optional<std::reference_wrapper<std::ostream>> out) = 0;
        virtual void load(std::ifstream& file) = 0;
        virtual void save(std::ofstream& file) = 0;
        virtual ~IEngine() = default;
    };

    template <typename p_t, typename vel_t, typename vflow_t, int valN, int valM>
    class Engine : public IEngine {
    private:
        int N = valN;
        int M = valM;

        Array<char, valN, valM> field{};
        VectorField<vel_t, valN, valM> velocity = {};
        VectorField<vflow_t, valN, valM> velocity_flow = {};
        Array<p_t, valN, valM> p{}, old_p{};
        Array<int64_t, valM, valM> dirs{};
        Array<int, valM, valM> last_use{};
        int UT = 0;
        p_t rho[256];

        std::tuple<vflow_t, bool, pair<int, int>> propagate_flow(int x, int y, vflow_t lim) {
            last_use[x][y] = UT - 1;
            vflow_t ret{};
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;

                // Added checks to the original algorithm
                if (nx < 0 || nx >= N || ny < 0 || ny >= M) continue;

                if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
                    vel_t cap = velocity.get(x, y, dx, dy);
                    vflow_t flow = velocity_flow.get(x, y, dx, dy);
                    if (fabs(flow - vflow_t(cap)) <= 0.0001) {
                        continue;
                    }
                    // assert(v >= velocity_flow.get(x, y, dx, dy));
                    vflow_t vp = std::min(lim, vflow_t(cap) - flow);
                    if (last_use[nx][ny] == UT - 1) {
                        velocity_flow.add(x, y, dx, dy, vp);
                        last_use[x][y] = UT;
                        // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
                        return {vp, true, {nx, ny}};
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
            return {ret, false, {0, 0}};
        }

        void propagate_stop(int x, int y, bool force = false) {
            if (!force) {
                for (auto [dx, dy] : deltas) {
                    int nx = x + dx, ny = y + dy;

                    // Added checks to the original algorithm
                    if (nx < 0 || nx >= N || ny < 0 || ny >= M) continue;

                    if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > 0ll) {
                        return;
                    }
                }
            }
            last_use[x][y] = UT;
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;

                // Added checks to the original algorithm
                if (nx < 0 || nx >= N || ny < 0 || ny >= M) continue;

                if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > 0ll) {
                    continue;
                }
                propagate_stop(nx, ny);
            }
        }

        vel_t move_prob(int x, int y) {
            vel_t sum{};
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;

                // Added checks to the original algorithm
                if (nx < 0 || nx >= N || ny < 0 || ny >= M) continue;

                if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                    continue;
                }
                vel_t v = velocity.get(x, y, dx, dy);
                if (v < 0ll) {
                    continue;
                }
                sum += v;
            }
            return sum;
        }

        void swap(int x1, int y1, int x2, int y2) {
            std::swap(field[x1][y1], field[x2][y2]);
            std::swap(p[x1][y1], p[x2][y2]);
            std::swap(velocity.v[x1][y1], velocity.v[x2][y2]);
        }

        bool propagate_move(int x, int y, bool is_first) {
            last_use[x][y] = UT - is_first;
            bool ret = false;
            int nx = -1, ny = -1;
            do {
                std::array<vel_t, deltas.size()> tres;
                vel_t sum{};
                for (size_t i = 0; i < deltas.size(); ++i) {
                    auto [dx, dy] = deltas[i];
                    int fx = x + dx, fy = y + dy;

                    // Added checks to the original algorithm
                    if (fx < 0 || fx >= N || fy < 0 || fy >= M) continue;

                    if (field[fx][fy] == '#' || last_use[fx][fy] == UT) {
                        tres[i] = sum;
                        continue;
                    }
                    vel_t v = velocity.get(x, y, dx, dy);
                    if (v < 0ll) {
                        tres[i] = sum;
                        continue;
                    }
                    sum += v;
                    tres[i] = sum;
                }

                if (sum == 0ll) {
                    break;
                }

                vel_t randNum = random01<vel_t>() * sum;
                size_t d = std::ranges::upper_bound(tres, randNum) - tres.begin();

                auto [dx, dy] = deltas[d];
                nx = x + dx;
                ny = y + dy;
                assert(velocity.get(x, y, dx, dy) > 0ll && field[nx][ny] != '#' && last_use[nx][ny] < UT);

                ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
            } while (!ret);
            last_use[x][y] = UT;
            for (auto [dx, dy] : deltas) {
                int fx = x + dx, fy = y + dy;

                // Added checks to the original algorithm
                if (fx < 0 || fx >= N || fy < 0 || fy >= M) continue;

                if (field[fx][fy] != '#' && last_use[fx][fy] < UT - 1 && velocity.get(x, y, dx, dy) < 0ll) {
                    propagate_stop(nx, ny);
                }
            }
            if (ret && !is_first) {
                swap(x, y, nx, ny);
            }
            return ret;
        }

        void init() {
            velocity_flow.v.init(N, M);
            dirs.init(N, M);
            old_p.init(N, M);

            rho[' '] = 0.01;
            rho['.'] = 1000ll;
            for (int x = 0; x < N; ++x) {
                for (int y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    for (auto [dx, dy]: deltas) {
                        dirs[x][y] += (field[x + dx][y + dy] != '#');
                    }
                }
            }
        }

    public:
        Engine() = default;

        void next(std::optional<std::reference_wrapper<std::ostream>> out) override {
            p_t total_delta_p = 0ll;

            // Apply external forces
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    if (field[x + 1][y] != '#')
                        velocity.add(x, y, 1, 0, g<vel_t>());
                }
            }

            // Apply forces from p
            old_p = p;
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    for (auto [dx, dy] : deltas) {
                        int nx = x + dx, ny = y + dy;
                        if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
                            p_t delta_p = old_p[x][y] - old_p[nx][ny];
                            p_t force = delta_p;
                            vel_t &contr = velocity.get(nx, ny, -dx, -dy);
                            if (p_t(contr) * rho[(int) field[nx][ny]] >= force) {
                                contr -= vel_t(force / rho[(int) field[nx][ny]]);
                                continue;
                            }
                            force -= p_t(contr) * rho[(int) field[nx][ny]];
                            contr = 0ll;
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
                prop = false;
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] != '#' && last_use[x][y] != UT) {
                            auto [t, local_prop, _] = propagate_flow(x, y, 1ll);
                            if (t > 0ll) {
                                prop = true;
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
                        vel_t old_v = velocity.get(x, y, dx, dy);
                        vflow_t new_v = velocity_flow.get(x, y, dx, dy);
                        if (old_v > 0ll) {
                            assert(vel_t(new_v) <= old_v);

                            velocity.get(x, y, dx, dy) = vel_t(new_v);
                            auto force = p_t(old_v - vel_t(new_v)) * rho[(int) field[x][y]];
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
                        if (random01<vel_t>() < move_prob(x, y)) {
                            prop = true;
                            propagate_move(x, y, true);
                        } else {
                            propagate_stop(x, y, true);
                        }
                    }
                }
            }

            if (prop && out) {
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < M; ++j) {
                        out->get() << field[i][j];
                    }
                    out->get() << std::endl;
                }
                std::cout.flush();
            }
        }

        void load(std::ifstream& file) override {
            auto loadArray = [&]<typename T, int N, int M>(Array<T, N, M>& arr, int n, int m) {
                arr.init(n, m);
                for (int i = 0; i < arr.N; ++i) {
                    for (int j = 0; j < arr.M; ++j) {
                        if constexpr (std::is_same_v<T, char>) {
                            char tmp;
                            do {
                                file.get(tmp);
                            } while (tmp == '\n');
//                            file >> tmp;
                            arr[i][j] = static_cast<char>(tmp);
                        } else {
                            double tmp = 0;
                            file >> tmp;
                            arr[i][j] = tmp;
                        }
                    }
                }
            };

            auto loadArrayOfArrays = [&]<typename T, int N, int M>(Array<std::array<T, 4>, N, M>& arr, int n, int m){
                arr.init(n, m);
                for (int i = 0; i < arr.N; ++i) {
                    for (int j = 0; j < arr.M; ++j) {
                        double tmp[4] = {0};
                        file >> tmp[0] >> tmp[1] >> tmp[2] >> tmp[3];
                        arr[i][j] = {T(tmp[0]), T(tmp[1]), T(tmp[2]), T(tmp[3])};
                    }
                }
            };


            if (!file.is_open()) {
                throw std::invalid_argument("file is not open");
            }
            file >> N >> M >> UT;
            loadArray(field, N, M);
            loadArray(last_use, N, M);
            loadArray(p, N, M);
            loadArrayOfArrays(velocity.v, N, M);
            init();
        }

        void save(std::ofstream &file) override {
            auto saveArray = [&]<typename T, int N, int M>(Array<T, N, M>& arr) {
                for (int i = 0; i < arr.N; ++i) {
                    for (int j = 0; j < arr.M; ++j) {
                        if constexpr (std::is_same_v<T, char>) {
                            file << arr[i][j];
                        } else {
                            file << arr[i][j] << " ";
                        }
                    }
                    file << std::endl;
                }
            };

            auto saveArrayOfArrays = [&]<typename T, int N, int M>(Array<std::array<T, 4>, N, M>& arr){
                for (int i = 0; i < arr.N; ++i) {
                    for (int j = 0; j < arr.M; ++j) {
                        file << arr[i][j][0] << " " << arr[i][j][1] << " " << arr[i][j][2] << " " << arr[i][j][3];
                    }
                    file << std::endl;
                }
            };

            if (!file.is_open()) {
                throw std::invalid_argument("file is not open");
            }
            file << N << " " << M << " " << UT << std::endl;
            saveArray(field);
            saveArray(last_use);
            saveArray(p);
            saveArrayOfArrays(velocity.v);
        }
    };
}

#endif //FLUID_PROJECT_FLUIDENGINE_H
