#include <fdaPDE/core/fdaPDE/optimization.h>

#include <unsupported/Eigen/SparseExtra>

#include <filesystem>
using namespace fdapde;

int main(){
  
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>; 

    ScalarField<1> F1;
    F1 = [](Eigen::Matrix<double,1,1> x) -> double {
        return x[0]*(x[0]-1.0);
    };

    ScalarField<2> F2;
    F2= [](Eigen::Matrix<double,2,1> x) -> double {
        return (x[0]-0.25)*(x[0]-0.25)+(x[1]+0.25)*(x[1]+0.25);
    };

    ScalarField<1> F1b;
    F1b = [](Eigen::Matrix<double,1,1> x) -> double {
        return (x[0]-3.0)*(x[0]-3.0);
    };

    ScalarField<2> F2b;
    F2b = [](Eigen::Matrix<double,2,1> x) -> double {
        return (x[0]-2)*(x[0]-2)+(x[1]+2)*(x[1]+2);
    };
    
    vector_t x0_1 = 0.01*vector_t::Ones(1);
    vector_t x0_2 = 0.01*vector_t::Ones(2);
    vector_t x0_1b = 0.01*vector_t::Ones(1);
    vector_t x0_2b = 0.01*vector_t::Ones(2);

    vector_t l_1 = -1*vector_t::Ones(1);
    vector_t u_1 = 1*vector_t::Ones(1);
    vector_t l_2 = -2*vector_t::Ones(2);
    vector_t u_2 = 2*vector_t::Ones(2);
    
    std::cout << "Testing of the LBFGSB method implementation" << std::endl;

    /////////////////////////////////////////
    // UNBOUNDED OPTIMIZATION IN 1D AND 2D //
    /////////////////////////////////////////

    std::cout << "\t LBFGSB (1)" << std::endl;
    LBFGSB<1> lbfgsb1;
    lbfgsb1.optimize(F1,x0_1, BacktrackingLineSearch());
    std::cout << "x opt: " << lbfgsb1.optimum()[0]  << std::endl;
    std::cout << "value opt: " << lbfgsb1.value() << std::endl;
    std::cout << "#iters: " << lbfgsb1.n_iter() << std::endl;
    // Correct result : x opt: 0.5, value opt: -0.25

    std::cout << "\t LBFGSB (2)" << std::endl;
    LBFGSB<2> lbfgsb2;
    lbfgsb2.optimize(F2, x0_2, BacktrackingLineSearch());
    std::cout << "x opt: [" << lbfgsb2.optimum()[0] << ", " << lbfgsb2.optimum()[1] << " ]"  << std::endl;
    std::cout << "value opt: " << lbfgsb2.value() << std::endl;
    std::cout << "#iters: " << lbfgsb2.n_iter() << std::endl;
    // Correct result: x opt: (0.25 , -0.25), value opt: 0
    
    ///////////////////////////////////////
    // BOUNDED OPTIMIZATION IN 1D AND 2D //
    ///////////////////////////////////////

    std::cout << "\t LBFGSB (1b)" << std::endl;
    LBFGSB<1> lbfgsb1b;
    lbfgsb1b.optimize(F1b,x0_1b, l_1, u_1, BacktrackingLineSearch());
    std::cout << "x opt: " << lbfgsb1b.optimum()[0]  << std::endl;
    std::cout << "value opt: " << lbfgsb1b.value() << std::endl;
    std::cout << "#iters: " << lbfgsb1b.n_iter() << std::endl;
    // Correct result : x opt: 1, value opt: 4

    std::cout << "\t LBFGSB (2b)" << std::endl;
    LBFGSB<2> lbfgsb2b;
    lbfgsb2b.optimize(F2b, x0_2b, l_2, u_2, BacktrackingLineSearch());
    std::cout << "x opt: [" << lbfgsb2b.optimum()[0] << ", " << lbfgsb2b.optimum()[1] << " ]"  << std::endl;
    std::cout << "value opt: " << lbfgsb2b.value() << std::endl;
    std::cout << "#iters: " << lbfgsb2b.n_iter() << std::endl;
    // Correct result: x opt: (2, -2), value opt: 0

    ///////////////////////////////////////
    // SOME CHALLENGING CASES (BOUNDED)  //
    ///////////////////////////////////////

    ScalarField<1> F1c;
    F1c = [](Eigen::Matrix<double,1,1> x) -> double {
        return (x[0]-0.001)*(x[0]-0.001) + std::exp(-1000*(x[0]-0.001)*(x[0]-0.001));
    };

    ScalarField<2> F2c; // Rosenbrock function
    F2c= [](Eigen::Matrix<double,2,1> x) -> double {
        return (1-x[0])*(1-x[0]) + 100*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0]);
    };

    vector_t x0_1c = 0.5*vector_t::Ones(1);
    vector_t x0_2c = -1.5*vector_t::Ones(2);

    vector_t l_1c = -1*vector_t::Ones(1);
    vector_t u_1c = 1*vector_t::Ones(1);
    vector_t l_2c = -2*vector_t::Ones(2);
    vector_t u_2c = 2*vector_t::Ones(2);

    std::cout << "\t LBFGSB (1c)" << std::endl;
    LBFGSB<1> lbfgsb1c;
    lbfgsb1c.optimize(F1c,x0_1c, l_1c, u_1c, BacktrackingLineSearch());
    std::cout << "x opt: " << lbfgsb1c.optimum()[0]  << std::endl;
    std::cout << "value opt: " << lbfgsb1c.value() << std::endl;
    std::cout << "#iters: " << lbfgsb1c.n_iter() << std::endl;
    // Correct result : x opt: 0.084113, value opt: 0.00791

    std::cout << "\t LBFGSB (2c)" << std::endl;
    LBFGSB<2> lbfgsb2c;
    lbfgsb2c.optimize(F2c, x0_2c, l_2c, u_2c, BacktrackingLineSearch());
    std::cout << "x opt: [" << lbfgsb2c.optimum()[0] << ", " << lbfgsb2c.optimum()[1] << " ]"  << std::endl;
    std::cout << "value opt: " << lbfgsb2c.value() << std::endl;
    std::cout << "#iters: " << lbfgsb2c.n_iter() << std::endl;
    // Correct result: x opt: (1, -1), value opt: 0
    
    return 0;
}
