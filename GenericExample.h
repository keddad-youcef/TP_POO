#include <iostream>
#include "ResNorm.h"
#include "AlienMock.h"
#include <memory>
// TD4 //


enum class SolverType
{
    Hypre,
    PETSc
};


class GenericExample {
    public:
    int run(SolverType m);

    template <typename C>
    void info();
    void run_parallel_thread();
    
};


template <typename C>
void GenericExample::info() {
    C::info();
}



class UniqueAPI {
public : 
  virtual void solve() = 0;
  virtual std::unique_ptr<Alien::ILinearAlgebra> createAlgebra() = 0;
  virtual std::unique_ptr<Alien::ILinearSolver> createSolver() = 0;
  virtual ~UniqueAPI() = default;
};

class HypreAPI : public UniqueAPI {
public:
    void solve() override {std::cout << "HypreAPI solver" << std::endl;}
    std::unique_ptr<Alien::ILinearAlgebra> createAlgebra() override;
    std::unique_ptr<Alien::ILinearSolver> createSolver() override;

    static void info();
};

class PETScAPI : public UniqueAPI {
public:
    void solve() override {std::cout << "PETScAPI solver" << std::endl;}
    std::unique_ptr<Alien::ILinearAlgebra> createAlgebra() override;
    std::unique_ptr<Alien::ILinearSolver> createSolver() override;

    static void info();
};


class SolverFabric {
    public:
    static auto create(SolverType st);
};
