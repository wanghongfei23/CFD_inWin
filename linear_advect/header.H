#include <algorithm>
#include <array>
#include <cmath>
#include <format> // Include std::format
#include <fstream>
#include <iomanip> // For std::setw
#include <iostream>
#include <omp.h>
#include <quadmath.h>
#include <stdexcept>
#include <vector>

typedef __float128 real_t;

class GhostArray {
private:
    std::vector<real_t> data;
    int ghost_size; // declared first
    int domain_size; // declared second

public:
    // 原有构造函数
    GhostArray(int domain_points, int ghost_points)
        : data(domain_points + 2 * ghost_points) // 必须先初始化data
        , ghost_size(ghost_points) // 然后ghost_size
        , domain_size(domain_points) // 最后domain_size
    {
    }

    // 修正后的拷贝构造函数
    GhostArray(const GhostArray& other)
        : data(other.data) // 保持与声明顺序一致
        , ghost_size(other.ghost_size)
        , domain_size(other.domain_size)
    {
    }

    // 新增获取ghost_size的公共方法
    int get_ghost_size() const { return ghost_size; }
    int get_domain_size() const { return domain_size; }

    // Total size including ghosts
    int total_size() const { return data.size(); }

    // Access operators with bounds checking
    real_t& operator[](int i)
    {
        if (i < -ghost_size || i >= domain_size + ghost_size) {
            throw std::out_of_range("Index out of ghost array bounds");
        }
        return data[i + ghost_size];
    }

    const real_t& operator[](int i) const
    {
        if (i < -ghost_size || i >= domain_size + ghost_size) {
            throw std::out_of_range("Index out of ghost array bounds");
        }
        return data[i + ghost_size];
    }

    // Direct access to underlying data
    real_t* raw_data() { return data.data(); }
    const real_t* raw_data() const { return data.data(); }

    // Fill ghost points for periodic boundaries
    void fill_periodic()
    {
        // Left ghost points
        for (int i = -ghost_size; i < 0; ++i) {
            data[i + ghost_size] = data[domain_size + i + ghost_size];
        }
        // Right ghost points
        for (int i = domain_size; i < domain_size + ghost_size; ++i) {
            data[i + ghost_size] = data[i - domain_size + ghost_size];
        }
    }

    // Get range of domain indices (excluding ghosts)
    struct Range {
        struct Iterator {
            int current;
            Iterator(int pos)
                : current(pos)
            {
            }

            int operator*() const { return current; }
            Iterator& operator++()
            {
                ++current;
                return *this;
            }
            Iterator operator++(int)
            {
                Iterator temp = *this;
                ++(*this);
                return temp;
            }
            Iterator& operator--()
            {
                --current;
                return *this;
            }
            Iterator operator--(int)
            {
                Iterator temp = *this;
                --(*this);
                return temp;
            }
            Iterator& operator+=(int n)
            {
                current += n;
                return *this;
            }
            Iterator operator+(int n) const
            {
                Iterator temp = *this;
                temp += n;
                return temp;
            }
            Iterator& operator-=(int n)
            {
                current -= n;
                return *this;
            }
            Iterator operator-(int n) const
            {
                Iterator temp = *this;
                temp -= n;
                return temp;
            }
            int operator-(const Iterator& other) const
            {
                return current - other.current;
            }
            bool operator==(const Iterator& other) const { return current == other.current; }
            bool operator!=(const Iterator& other) const { return !(*this == other); }
            bool operator<(const Iterator& other) const { return current < other.current; }
            bool operator>(const Iterator& other) const { return current > other.current; }
            bool operator<=(const Iterator& other) const { return !(*this > other); }
            bool operator>=(const Iterator& other) const { return !(*this < other); }
        };

        Iterator begin() const { return Iterator(start); }
        Iterator end() const { return Iterator(start + length); }

        int start;
        int length;
    };

    Range domain() const { return { 0, domain_size }; }
    Range full_range() const { return { -ghost_size, domain_size + 2 * ghost_size }; }
};

inline real_t Teno5_CongZ(const std::array<real_t, 5>& q)
{
    // real_t eps = 1.0e-40Q;
    std::array<real_t, 3> beta = {
        1.0Q / 1.0Q * powq(1.0Q * q[0] - 2.0Q * q[1] + 1.0Q * q[2], 2) + 1.0Q / 4.0Q * powq(1.0Q * q[0] - 4.0Q * q[1] + 3.0Q * q[2], 2),
        1.0Q / 1.0Q * powq(1.0Q * q[1] - 2.0Q * q[2] + 1.0Q * q[3], 2) + 1.0Q / 4.0Q * powq(1.0Q * q[1] + 0.0Q * q[2] - 1.0Q * q[3], 2),
        1.0Q / 1.0Q * powq(1.0Q * q[2] - 2.0Q * q[3] + 1.0Q * q[4], 2) + 1.0Q / 4.0Q * powq(3.0Q * q[2] - 4.0Q * q[3] + 1.0Q * q[4], 2)
    };
    unsigned short minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin();
    constexpr real_t CT = std::pow(1.5 * 1.0e-5, 1.0 / 6.0);
    constexpr real_t CT_1 = 1.0 - CT;
    real_t tau = std::abs(beta[2] - beta[0]);
    real_t rr = CT * tau - CT_1 * beta[minBeta];
    real_t ll = tau * beta[minBeta];
    unsigned short flag = 0;
    if (ll < rr * beta[0])
        flag += 1;
    if (ll < rr * beta[1])
        flag += 2;
    if (ll < rr * beta[2])
        flag += 4;
    switch (flag) {
    // switch (0) {
    case 0:
        /* 1,1,1 */
        return 3.0Q / 128.0Q * q[0] - 5.0Q / 32.0Q * q[1] + 45.0Q / 64.0Q * q[2] + 15.0Q / 32.0Q * q[3] - 5.0Q / 128.0Q * q[4];
    case 1:
        /* 0,1,1 */
        return -1.0Q / 16.0Q * q[1] + 9.0Q / 16.0Q * q[2] + 9.0Q / 16.0Q * q[3] - 1.0Q / 16.0Q * q[4];
    case 2:
        /* 1,0,1 */
        return 3.0Q / 8.0Q * q[2] + 3.0Q / 4.0Q * q[3] - 1.0Q / 8.0Q * q[4];
    case 3:
        /* 0,0,1 */
        return 3.0Q / 8.0Q * q[2] + 3.0Q / 4.0Q * q[3] - 1.0Q / 8.0Q * q[4];
    case 4:
        /* 1,1,0 */
        return 1.0Q / 16.0Q * q[0] - 5.0Q / 16.0Q * q[1] + 15.0Q / 16.0Q * q[2] + 5.0Q / 16.0Q * q[3];
    case 5:
        /* 0,1,0 */
        return -1.0Q / 8.0Q * q[1] + 3.0Q / 4.0Q * q[2] + 3.0Q / 8.0Q * q[3];
    case 6:
        /* 1,0,0 */
        return 3.0Q / 8.0Q * q[0] - 5.0Q / 4.0Q * q[1] + 15.0Q / 8.0Q * q[2];
    default:
        /* 0,0,0 */
        return q[2];
    }
}

real_t exact_initial_solution(real_t x, real_t /*t*/)
{
    // Example: exact solution for a periodic sine function
    // u(x, t) = sin(2 * pi * (x - ct))
    // real_t c = 1.0Q; // wave speed
    // real_t x0 = x - c * t;
    // return sinq(2.0Q * M_PIq * x0) + 1.0Q; // Use M_PIq for quad-precision pi

    // Exact solution for a Gaussian function
    // u(x, t) = e^(-300 * (x - x_c)^2) where x_c = 0.5
    // real_t x_c = 0.5Q;
    // return expq(-300.0Q * (x - x_c) * (x - x_c));
    return x * x * expq(0.75Q * (x - 1.0Q));
    // return x;
}
real_t exact_deriv(real_t x, real_t /*t*/)
{
    // real_t x_c = 0.5Q;
    // return 600.0Q * (x - x_c) * expq(-300.0Q * (x - x_c) * (x - x_c));
    // return -1.Q;
    return -expq(0.75Q * (x - 1.0Q)) * (2.0Q * x + 0.75Q * x * x);
}

real_t exact_solution(real_t x, real_t /*t*/)
{
    // return exact_initial_solution(x, 0.0Q);
    return exact_deriv(x, 0.0Q);
}

void print_real(real_t value)
{
    // Buffer for quadmath conversion
    char buf[128];

    // Convert quad precision to string with 10 decimal places
    quadmath_snprintf(buf, sizeof(buf), "%.10Qf", value);

    // Print to standard output
    std::cout << buf;
}
std::ostream& operator<<(std::ostream& os, const real_t& value)
{
    char buffer[128];
    quadmath_snprintf(buffer, sizeof(buffer), "%.30Qg", value);
    os << buffer;
    return os;
}

static bool is_periodic = true;

// 计算右端项 (WCNS5空间离散)
void compute_rhs(GhostArray& rhs, GhostArray& u, GhostArray& flux, real_t dx)
{
    const real_t c = 1.0Q; // 波速
    if (is_periodic)
        u.fill_periodic();

    // 1. 一阶插值计算界面通量 (i+1/2位置)
#pragma omp parallel for
    for (int i : flux.full_range()) {
        // 提取左值和右值
        // real_t u_left = u[i - 1];
        // real_t u_right = u[i];
        std::array<real_t, 5> q_left = { u[i - 3], u[i - 2], u[i - 1], u[i], u[i + 1] };
        std::array<real_t, 5> q_right = { u[i + 2], u[i + 1], u[i], u[i - 1], u[i - 2] };

        // 使用高阶插值得到左值和右值
        real_t u_left = Teno5_CongZ(q_left);
        real_t u_right = Teno5_CongZ(q_right);

        // 一阶插值:取左右平均值
        // if (c > 0) {
        //     // c > 0, 使用前向差分
        //     flux[i] = c * u_left;
        // } else {
        //     // c <= 0, 使用后向差分
        //     flux[i] = c * u_right;
        // }
        real_t F_left = c * u_left;
        real_t F_right = c * u_right;
        flux[i] = (F_left + F_right) / 2.0Q - (c / 2.0Q) * (u_right - u_left);
    }
    if (is_periodic)
        flux.fill_periodic();
    std::array<real_t, 3> w = { 3.0Q / 2.0Q, -3.0Q / 10.0Q, 1.0Q / 30.0Q };
// std::array<real_t, 3> w = { 4.0Q / 3.0Q, -1.0Q / 6.0Q, 0.0Q };
// std::array<real_t, 3> w = { 1.0Q, 0.0Q, 0.0Q };
// 2. 二阶中心差分计算右端项
#pragma omp parallel for
    for (int i : rhs.domain()) {
        // 二阶差分: (flux[i+1] - flux[i]) / dx
        // rhs[i] = -1.0Q / dx * (flux[i + 1] - flux[i]);
        rhs[i] = -w[0] / dx * (flux[i + 1] - flux[i])
            - w[1] / dx * (u[i + 1] - u[i - 1])
            - w[2] / dx * (flux[i + 2] - flux[i - 1]);
    }
}

void evolve(GhostArray& u, GhostArray& u_prev, GhostArray& flux, GhostArray& rhs,
    real_t dt, real_t dx) // 添加dx参数
{
    // 阶段1
    compute_rhs(rhs, u_prev, flux, dx); // 现在dx已定义
#pragma omp parallel for
    for (int i : u.domain()) {
        u[i] = u_prev[i] + dt * rhs[i];
    }

    // 阶段2
    GhostArray temp(u);
    compute_rhs(rhs, u, flux, dx); // 传递dx
#pragma omp parallel for
    for (int i : u.domain()) {
        temp[i] = 0.75Q * u_prev[i] + 0.25Q * (u[i] + dt * rhs[i]);
    }

    // 阶段3
    compute_rhs(rhs, temp, flux, dx); // 传递dx
#pragma omp parallel for
    for (int i : u.domain()) {
        u[i] = (1.0Q / 3.0Q) * u_prev[i] + (2.0Q / 3.0Q) * (temp[i] + dt * rhs[i]);
    }
}

// Function to calculate L1 norm of error
inline real_t L1_norm_error(const GhostArray& u, const GhostArray& x, real_t t)
{
    real_t error = 0.0Q;
    for (int i : u.domain()) {
        real_t xi = x[i];
        error += fabsq(u[i] - exact_solution(xi, t));
    }
    return error / u.get_domain_size();
}

// Function to calculate L2 norm of error
inline real_t L2_norm_error(const GhostArray& u, const GhostArray& x, real_t t)
{
    real_t error = 0.0Q;
    for (int i : u.domain()) {
        real_t xi = x[i];
        error += powq(u[i] - exact_solution(xi, t), 2);
    }
    return sqrtq(error / u.get_domain_size());
}

// Function to calculate L-infinity norm of error
inline real_t Linf_norm_error(const GhostArray& u, const GhostArray& x, real_t t)
{
    real_t max_error = 0.0Q;
    for (int i : u.domain()) {
        real_t xi = x[i];
        real_t error = fabsq(u[i] - exact_solution(xi, t));
        if (error > max_error) {
            max_error = error;
        }
    }
    return max_error;
}