#include <fdaPDE/core/fdaPDE/optimization.h>

#include <unsupported/Eigen/SparseExtra>

#include <filesystem>
using namespace fdapde;

int main(){
    using vector_t = Eigen::Matrix<double, Dynamic, 1>; 
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>; 

    ScalarField<1> F;
    F = [](Eigen::Matrix<double,1,1> x) -> double {
        return x[0]*(x[0]-1.0);
    };

    ScalarField<1> F2;
    F2 = [](Eigen::Matrix<double,1,1> x) -> double {
        return (x[0]-1.)*(x[0]-1.0);
    };

    int n_x = 100;
    vector_t x_grid(n_x);
    for (int i = 0; i < n_x; ++i) {
      x_grid(i, 0) = static_cast<double>(i+1)/n_x; 
    }
    
    std::cout << "\t Grid Search " << std::endl;
    GridSearch<1> optimizer;
    optimizer.optimize(F, x_grid);
    std::cout << "x opt: " << optimizer.optimum()[0]  << std::endl;
    std::cout << "value opt: " << optimizer.value() << std::endl;

    for(int i=0; i < optimizer.values().size(); ++i){
      std::cout << optimizer.values()[i] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "\t Grid Search " << std::endl;
    optimizer.optimize(F2, x_grid);
    std::cout << "x opt: " << optimizer.optimum()[0]  << std::endl;
    std::cout << "value opt: " << optimizer.value() << std::endl;

    for(int i=0; i < optimizer.values().size(); ++i){
      std::cout << optimizer.values()[i] << " ";
    }
    std::cout << std::endl;
    
    
    vector_t x0 = 0.01*vector_t::Ones(1);
    std::cout << "\t LBFGS " << std::endl;
    LBFGS<1> lbfgs;
    lbfgs.optimize(F,x0, BacktrackingLineSearch());
    std::cout << "x opt: " << lbfgs.optimum()[0]  << std::endl;
    std::cout << "value opt: " << lbfgs.value() << std::endl;
    std::cout << "#iters: " << lbfgs.n_iter() << std::endl;
    
    std::cout << "\t BFGS " << std::endl;
    BFGS<1> bfgs;
    bfgs.optimize(F,x0, BacktrackingLineSearch());
    std::cout << "x opt: " << bfgs.optimum()[0]  << std::endl;
    std::cout << "value opt: " << bfgs.value() << std::endl;
    std::cout << "#iters: " << bfgs.n_iter() << std::endl;
    
    std::cout << "\t Nelder-Mead " << std::endl;
    NelderMead<1> nelder_mead;
    nelder_mead.optimize(F,x0, BacktrackingLineSearch());
    std::cout << "x opt: " << nelder_mead.optimum()[0]  << std::endl;
    std::cout << "value opt: " << nelder_mead.value() << std::endl;
    std::cout << "#iters: " << nelder_mead.n_iter() << std::endl;
    

    return 0;
}
