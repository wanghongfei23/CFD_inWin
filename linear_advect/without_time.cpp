#include "header.H"

struct GridVariables {
    GhostArray x; // 网格坐标
    GhostArray u; // 解变量
    GhostArray u_prev; // 前一时间步的解
    GhostArray flux; // 通量存储在半节点上 (i+1/2位置)
    GhostArray rhs; // 右端项

    real_t dx; // 空间步长
    real_t time; // 当前时间

    // 正确的构造函数
    GridVariables(int N, int ghost, real_t length)
        : x(N, ghost)
        , // N个网格点
        u(N, ghost)
        , u_prev(N, ghost)
        , flux(N + 1, ghost - 3)
        , // N个区间有N+1个界面
        rhs(N, 0)
        , // 右端项与解变量同尺寸
        dx(length / N)
        , time(0.0Q)
    {
        // 初始化网格坐标 (可选)
        // for (int i : x.full_range()) {
        //     x[i] = i * dx;
        // }
        // x.fill_periodic();
    }
};

void write_to_tecplot_file(const GridVariables& grid, real_t current_time, const std::string& filename)
{
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        throw std::runtime_error("Failed to open file for writing.");
    }

    outfile << "TITLE = \"Simulation Data\"\n";
    outfile << "VARIABLES = \"X\", \"U\"\n";
    outfile << "ZONE T=\"t=" << current_time << "\", I=" << grid.x.domain().length << "\n";

    for (int i : grid.x.domain()) {
        outfile << grid.x[i] << " " << grid.u[i] << "\n";
    }

    outfile.close();
}
int main()
{
    omp_set_num_threads(10);

    is_periodic = false;
    // Grid density sequence
    std::vector<int> grid_densities = { 50, 100, 150, 200, 300, 400, 500 };

    // Parameters
    const int ghost = 5; // 每侧虚拟点数
    const real_t L = 1.0Q; // 域长度

    // Vectors to store errors for each grid density
    std::vector<real_t> L1_errors;
    std::vector<real_t> L2_errors;
    std::vector<real_t> Linf_errors;

    std::ofstream results_file("rhs_results.txt");
    if (!results_file.is_open()) {
        std::cerr << "Failed to open rhs_results.txt for writing." << std::endl;
        return 1;
    }

    // Write headers to the file
    results_file << "N\tL1_Error\tL2_Error\tLinf_Error\tOrder_of_Accuracy\n";

    for (int N : grid_densities) {
        // Create grid variables collection
        GridVariables grid(N, ghost, L);

        // Initialize grid
#pragma omp parallel for
        for (GhostArray::Range::Iterator i = grid.x.full_range().begin(); i != grid.x.full_range().end(); ++i) {
            grid.x[*i] = (*i + 0.5Q) * grid.dx; // 单元中心型网格坐标
            // 方波初始条件
            grid.u[*i] = exact_initial_solution(grid.x[*i], 0.0);
        }

        // Compute RHS
        compute_rhs(grid.rhs, grid.u, grid.flux, grid.dx);

        // Calculate errors
        real_t L1_error = L1_norm_error(grid.rhs, grid.x, 0.0Q);
        real_t L2_error = L2_norm_error(grid.rhs, grid.x, 0.0Q);
        real_t Linf_error = Linf_norm_error(grid.rhs, grid.x, 0.0Q);

        // Store errors
        L1_errors.push_back(L1_error);
        L2_errors.push_back(L2_error);
        Linf_errors.push_back(Linf_error);

        // Write errors to file with initial order of accuracy set to zero
        // results_file << N << "\t" << L1_error << "\t" << L2_error << "\t" << Linf_error << "\t0\n";

        // Output errors
        std::cout << "Grid density: " << N << "\n";
        std::cout << "L1 norm error: " << L1_error << "\n";
        std::cout << "L2 norm error: " << L2_error << "\n";
        std::cout << "L-infinity norm error: " << Linf_error << "\n\n";
    }

    // Calculate order of accuracy
    std::cout << "Order of accuracy:\n";
    results_file << std::format("{:d}\t{:.5g}\t0\t{:.5g}\t0\t{:.5g}\t0\n",
        grid_densities[0], L1_errors[0], L2_errors[0], Linf_errors[0]);
    for (size_t i = 1; i < grid_densities.size(); ++i) {
        real_t h1 = 1.0Q / grid_densities[i - 1];
        real_t h2 = 1.0Q / grid_densities[i];

        real_t L1_order = log2(L1_errors[i - 1] / L1_errors[i]) / log2(h1 / h2);
        real_t L2_order = log2(L2_errors[i - 1] / L2_errors[i]) / log2(h1 / h2);
        real_t Linf_order = log2(Linf_errors[i - 1] / Linf_errors[i]) / log2(h1 / h2);

        std::cout << "From N = " << grid_densities[i - 1] << " to N = " << grid_densities[i] << ":\n";
        std::cout << "L1 norm order: " << L1_order << "\n";
        std::cout << "L2 norm order: " << L2_order << "\n";
        std::cout << "L-infinity norm order: " << Linf_order << "\n\n";

        // Update the order of accuracy
        results_file << std::format("{:d}\t{:.5g}\t{:.5g}\t{:.5g}\t{:.5g}\t{:.5g}\t{:.5g}\n",
            grid_densities[i], L1_errors[i], L1_order,
            L2_errors[i], L2_order, Linf_errors[i], Linf_order);
    }
    results_file.close();
    return 0;
}