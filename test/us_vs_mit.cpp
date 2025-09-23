// CODE TO TEST THE PERFORMANCE OF OUR LBFGSB OPTIMIZER IMPLEMENTATION
// AGAINST MIT'S ONE (2D CASE)

// Make necessary includes
#include <LBFGSB.h>                                     // MIT version
#include <fdaPDE/core/fdaPDE/optimization.h>            // OUR version

#include <LBFGSpp/LineSearchBacktracking.h>
#include <fdaPDE/fdapde.h>
#include <unsupported/Eigen/SparseExtra>
#include <filesystem>
#include <chrono>

// useful type aliases
using Eigen::VectorXd;
using fdapde::ScalarField;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

int main() {
    // Define MIT minimizer
    LBFGSpp::LBFGSBParam<double> param;
    param.epsilon = 1e-6;
    param.max_iterations = 1e4;
    LBFGSpp::LBFGSBSolver<double> solver_mit(param);
    // Define OUR minimzer
    fdapde::LBFGSB<2> solver_our;

    // Define bounds vector (for constrained optimization)
    int min = -100, max = 100;
    VectorXd lb = VectorXd::Constant(2, min);
    VectorXd ub = VectorXd::Constant(2, max);

    // Define the function to minimize (a simple paraboloid in 2D in this case)
    ScalarField<2> fun;
    fun = [](Eigen::Matrix<double,2,1> x) -> double {
        return x[0]*x[0] + x[1]*x[1];
    };
    // MIT implementation requires to create this type of struct starting from the scalar field
    struct function_to_optimize {
        function_to_optimize(const ScalarField<2>& Field) : Field_(Field) { }

        double operator()(const vector_t& x, vector_t& grad) {
            double Fx = Field_(x);
            grad = Field_.gradient()(x);
            return Fx;
        }

        const ScalarField<2>& Field_;
    };

    function_to_optimize fun_mit(fun);

    // Perform multiple simulations
    int num_simulations = 5;
    for(int i=0; i < num_simulations; ++i) {
        // Random initial guess
        Eigen::VectorXd x0(2);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(min, max);

        x0(0) = dis(gen);
        x0(1) = dis(gen);

        // Now perform the minimization
        std::cout << "________________________________________________________________" << std::endl;
        std::cout << "TEST NUMBER : " << i+1 << std::endl;
        std::cout << "Initial guess = (" << x0(0) << ", " << x0(1) << ")" << std::endl;
        // MIT version
        std::cout << "________________________________________________________________" << std::endl;
        double fx;                      // stores f(x) at the minimum
        Eigen::VectorXd x = x0;         // necessary since the solver overwrites the initial guess
        // Start timer (MIT)
        auto start_mit = std::chrono::high_resolution_clock::now();                                     
        int niter = solver_mit.minimize(fun_mit, x, fx, lb, ub);
        // Stop timer (MIT)
        auto end_mit = std::chrono::high_resolution_clock::now();
        // Convert timer to nanoseconds (MIT)
        auto elapsed_ms_mit = std::chrono::duration_cast<std::chrono::nanoseconds>(end_mit - start_mit);
        // Output (MIT)
        std::cout << "MIT version" << std::endl;
        std::cout << niter << " iterations" << std::endl;
        std::cout << "x = (" << x(0) << ", " << x(1) << ")" << std::endl;
        std::cout << "f(x) = " << fx << std::endl;
        std::cout << "Elapsed time = " << elapsed_ms_mit << " nanoseconds" << std::endl; 
        // OUR version
        std::cout << "________________________________________________________________" << std::endl;
        // Start timer (OUR)
        auto start_our = std::chrono::high_resolution_clock::now();
        solver_our.optimize(fun, x0, fdapde::BacktrackingLineSearch());
        // Stop timer (OUR)
        auto end_our = std::chrono::high_resolution_clock::now();
        // Convert timer to nanoseconds (OUR)
        auto elapsed_ms_our = std::chrono::duration_cast<std::chrono::nanoseconds>(end_our - start_our);
        // Output (OUR)
        std::cout << "OUR version : " << std::endl;
        std::cout << solver_our.n_iter() << " iterations" << std::endl;
        std::cout << "x = (" << solver_our.optimum()[0] << ", " << solver_our.optimum()[1] << ")" << std::endl;
        std::cout << "f(x) = " << solver_our.value() << std::endl;
        std::cout << "Elapsed time = " << elapsed_ms_our << " nanoseconds" << std::endl; 
        std::cout << "________________________________________________________________" << std::endl;
    }
    // TO IMPLEMENT -> save everything to a csv file in order to make graphs easily
    return 0;
}