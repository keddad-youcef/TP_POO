#include <iostream>
#include <memory>
#include <thread>
#include "GenericExample.h"
#include "HypreExample.h"
#include "AlienMock.h"
#include "LocalLinearAlgebra.h"
#include "Timer.h"

// TD4 //

// SolverFabic create() method to create uinque_api smart pointer to either Hypre or PETSc.

auto SolverFabric::create(SolverType st) {
  
  // UniqueAPI* unique_api;
  // unique_api = new HypreAPI{};
  // delete unique_api;

  std::unique_ptr<UniqueAPI> unique_api;
  switch (st)
  {
  case SolverType::Hypre :
    unique_api = std::make_unique<HypreAPI>();
    // unique_api->solve();
    break;
  
  case SolverType::PETSc :
    unique_api = std::make_unique<PETScAPI>();
    // unique_api->solve();
    break;
  }
  return unique_api;
}



// createAlgebra() method: It creates a smart pointer to Hypre Algebra or PETSc Algebra. 

std::unique_ptr<Alien::ILinearAlgebra>  HypreAPI::createAlgebra() {
  std::unique_ptr<Alien::ILinearAlgebra> algebra_ptr = std::make_unique<Alien::Hypre::LinearAlgebra>();
  return algebra_ptr;
}

std::unique_ptr<Alien::ILinearAlgebra> PETScAPI::createAlgebra() {
  std::unique_ptr<Alien::ILinearAlgebra> algebra_ptr = std::make_unique<Alien::SimpleCSRLinearAlgebra>();
  return algebra_ptr;
}



// createSolver() method: It creates a smart pointer to Hypre Solver or PETSc Solver. 

std::unique_ptr<Alien::ILinearSolver> HypreAPI::createSolver() {
  std::unique_ptr<Alien::ILinearSolver> solver_ptr = std::make_unique<Alien::Hypre::LinearSolver>();
  return solver_ptr;
}


std::unique_ptr<Alien::ILinearSolver> PETScAPI::createSolver() {
  Alien::PETSc::Options options;
  options.numIterationsMax(100);
  options.stopCriteriaValue(1e-10);
  options.preconditioner(Alien::PETSc::OptionTypes::Jacobi);
  options.solver(Alien::PETSc::OptionTypes::BiCGstab /*CG*/);

  std::unique_ptr<Alien::ILinearSolver> solver_ptr = std::make_unique<Alien::PETSc::LinearSolver>();

  return solver_ptr;
}


// HypreAPI & PETScAPI info() methods.


void HypreAPI::info() {
  std::cout << "Library Name :  Hypre" << "\n\n";
}

void PETScAPI::info() {
  std::cout << "Library Name :  PETSc" << std::endl;
  std::cout << "preconditioner :  Jacobi" << std::endl;
  std::cout << "solverType :  BiCGstab" << std::endl;
  std::cout << "numIterationsMax =  100" << std::endl;
  std::cout << "stopCriteriaValue =  1e-10" << "\n\n";
}





int GenericExample::run(SolverType m) 
{
  
  auto* pm = Arccore::MessagePassing::Mpi::StandaloneMpiMessagePassingMng::create(MPI_COMM_WORLD);
  auto* tm = Arccore::arccoreCreateDefaultTraceMng();

  Alien::setTraceMng(tm);
  Alien::setVerbosityLevel(Alien::Verbosity::Debug);

  auto size = 100;

  tm->info() << "Example Alien :";
  tm->info() << "Use of scalar builder (RefSemantic API) for Laplacian problem";
  tm->info() << " => solving linear system Ax = b";
  tm->info() << " * problem size = " << size;
  tm->info() << " ";
  tm->info() << "Start example...";
  tm->info() << " ";

  Alien::Matrix A(size, size, pm);

  
  

  // Distributions calculée
  const auto& dist = A.distribution();
  int offset = dist.rowOffset();
  int lsize = dist.localRowSize();
  int gsize = dist.globalRowSize();

  tm->info() << "build matrix with direct matrix builder";
  {
    Alien::DirectMatrixBuilder builder(A, Alien::DirectMatrixOptions::eResetValues);
    builder.reserve(3); // Réservation de 3 coefficients par ligne
    builder.allocate(); // Allocation de l'espace mémoire réservé

    for (int irow = offset; irow < offset + lsize; ++irow) {
      builder(irow, irow) = 2.;
      if (irow - 1 >= 0)
        builder(irow, irow - 1) = -1.;
      if (irow + 1 < gsize)
        builder(irow, irow + 1) = -1.;
    }
  }

  tm->info() << "* xe = 1";

  Alien::Vector xe = Alien::ones(size, pm);

  tm->info() << "=> Vector Distribution : " << xe.distribution();

  tm->info() << "* b = A * xe";

  Alien::Vector b(size, pm);
  


  Alien::Vector x(size, pm);


  Alien::Vector r(size, pm);

  


  switch(m) {
  case SolverType::Hypre:
  {
    auto unique_api = SolverFabric::create(SolverType::Hypre);    // Create HypreAPI Pointer
    auto algebra = unique_api->createAlgebra();                   // Create Hypre algebra Pointer

    // Alien::Hypre::LinearAlgebra algebra;
    
    algebra->mult(A, xe, b);

    

    tm->info() << "* x = A^-1 b";

    //  auto options = Alien::Hypre::Options()
    //          .numIterationsMax(100)
    //          .stopCriteriaValue(1e-10)
    //          .preconditioner(Alien::Hypre::OptionTypes::AMGPC)
    //          .solver(Alien::Hypre::OptionTypes::GMRES);
    //
    //  auto solver = Alien::Hypre::LinearSolver (options);

    auto solver = unique_api->createSolver();                     // Create Hypre solver pointer

    // auto solver = Alien::Hypre::LinearSolver();

    solver->solve(A, b, x);

    tm->info() << "* r = Ax - b";

    

    {
      Alien::Vector tmp(size, pm);
      tm->info() << "t = Ax";
      algebra->mult(A, x, tmp);
      tm->info() << "r = t";
      algebra->copy(tmp, r);
      tm->info() << "r -= b";
      algebra->axpy(-1., b, r);
    }

    // auto 
    double norm = algebra->norm2(r);

    tm->info() << " => ||r|| = " << norm;

    tm->info() << "* r = || x - xe ||";

    {
      tm->info() << "r = x";
      algebra->copy(x, r);
      tm->info() << "r -= xe";
      algebra->axpy(-1., xe, r);
    }

    tm->info() << " => ||r|| = " << norm;

    tm->info() << " ";
    tm->info() << "... example finished !!!";


    std::cout << "Hypre Call" << "\n\n\n";
  }
  break;

  case SolverType::PETSc:
  {
    auto unique_api = SolverFabric::create(SolverType::PETSc);          // Create PETScAPI Pointer
    auto algebra = unique_api->createAlgebra();                         // Create PETSc algebra Pointer


    // Alien::SimpleCSRLinearAlgebra algebra;

    algebra->mult(A, xe, b);

   
  
    tm->info() << "* x = A^-1 b";

    Alien::PETSc::Options options;
    options.numIterationsMax(100);
    options.stopCriteriaValue(1e-10);
    options.preconditioner(Alien::PETSc::OptionTypes::Jacobi);
    options.solver(Alien::PETSc::OptionTypes::BiCGstab /*CG*/);
    //

    auto solver = unique_api->createSolver();                            // Create PETSc solver Pointer 

    // auto solver = Alien::PETSc::LinearSolver(options);
    //auto solver = Alien::PETSc::LinearSolver();

    
    solver->solve(A, b, x);

    tm->info() << "* r = Ax - b";


    {
      Alien::Vector tmp(size, pm);
      tm->info() << "t = Ax";
      algebra->mult(A, x, tmp);
      tm->info() << "r = t";
      algebra->copy(tmp, r);
      tm->info() << "r -= b";
      algebra->axpy(-1., b, r);
    }

    double norm = algebra->norm2(r);
    
    tm->info() << " => ||r|| = " << norm;



    tm->info() << "* r = || x - xe ||";

    {
      tm->info() << "r = x";
      algebra->copy(x, r);
      tm->info() << "r -= xe";
      algebra->axpy(-1., xe, r);
    }

    tm->info() << " => ||r|| = " << norm;

    tm->info() << " ";
    tm->info() << "... example finished !!!";

      std::cout << "Petsc Call" << "\n\n\n";
  }
    break;

  } 

  return 0;
};








//TD5



void GenericExample::run_parallel_thread() {
  std::cout << "Parallel Calcul" << "\n\n\n";

  
  
  std::thread first_thread([](){
    auto* pm = Arccore::MessagePassing::Mpi::StandaloneMpiMessagePassingMng::create(MPI_COMM_WORLD);
    auto* tm = Arccore::arccoreCreateDefaultTraceMng();

    Alien::setTraceMng(tm);
    Alien::setVerbosityLevel(Alien::Verbosity::Debug);

    auto size = 100;

    tm->info() << "Example Alien :";
    tm->info() << "Use of scalar builder (RefSemantic API) for Laplacian problem";
    tm->info() << " => solving linear system Ax = b";
    tm->info() << " * problem size = " << size;
    tm->info() << " ";
    tm->info() << "Start example...";
    tm->info() << " ";

    Alien::Matrix A(size, size, pm);
    

    // Distributions calculée
    const auto& dist = A.distribution();
    int offset = dist.rowOffset();
    int lsize = dist.localRowSize();
    int gsize = dist.globalRowSize();

    tm->info() << "build matrix with direct matrix builder";
    {
      Alien::DirectMatrixBuilder builder(A, Alien::DirectMatrixOptions::eResetValues);
      builder.reserve(3); // Réservation de 3 coefficients par ligne
      builder.allocate(); // Allocation de l'espace mémoire réservé

      for (int irow = offset; irow < offset + lsize; ++irow) {
        builder(irow, irow) = 2.;
        if (irow - 1 >= 0)
          builder(irow, irow - 1) = -1.;
        if (irow + 1 < gsize)
          builder(irow, irow + 1) = -1.;
      }
    }

    tm->info() << "* xe = 1";

    Alien::Vector xe = Alien::ones(size, pm);

    tm->info() << "=> Vector Distribution : " << xe.distribution();

    tm->info() << "* b = A * xe";

    Alien::Vector b(size, pm);

    Alien::Hypre::LinearAlgebra algebra;

    algebra.mult(A, xe, b);

    Alien::Vector x(size, pm);

    tm->info() << "* x = A^-1 b";

    //  auto options = Alien::Hypre::Options()
    //          .numIterationsMax(100)
    //          .stopCriteriaValue(1e-10)
    //          .preconditioner(Alien::Hypre::OptionTypes::AMGPC)
    //          .solver(Alien::Hypre::OptionTypes::GMRES);
    //
    //  auto solver = Alien::Hypre::LinearSolver (options);

    auto solver = Alien::Hypre::LinearSolver();

    solver.solve(A, b, x);

    tm->info() << "* r = Ax - b";

    Alien::Vector r(size, pm);

    {
      Alien::Vector tmp(size, pm);
      tm->info() << "t = Ax";
      algebra.mult(A, x, tmp);
      tm->info() << "r = t";
      algebra.copy(tmp, r);
      tm->info() << "r -= b";
      algebra.axpy(-1., b, r);
    }

    auto norm = algebra.norm2(r);

    tm->info() << " => ||r|| = " << norm;

    tm->info() << "* r = || x - xe ||";

    {
      tm->info() << "r = x";
      algebra.copy(x, r);
      tm->info() << "r -= xe";
      algebra.axpy(-1., xe, r);
    }

    tm->info() << " => ||r|| = " << norm;

    tm->info() << " ";
    tm->info() << "... example finished !!!";
    
  });


  std::thread second_thread([](){
    auto size = 100;
    LocalLinearAlgebra::Vector tmp_local(size);
    
    LocalLinearAlgebra::Vector r_local;

    LocalLinearAlgebra::Vector b_local(size);

    // dimensionnement et initialisation à 1 du vecteur
    LocalLinearAlgebra::Vector x_local(size, 1);
    
    // dimensionnement et initialisation à 0 de la matrice
    LocalLinearAlgebra::Matrix A_local(size, size, 0);

    LocalLinearAlgebra::mult(A_local, x_local, tmp_local); 
    
    LocalLinearAlgebra::mult(A_local, x_local, b_local);

    
    //
    LocalLinearAlgebra::copy(tmp_local, r_local);

    //
    LocalLinearAlgebra::axpy(-1., b_local, r_local);



    auto norm_local = LocalLinearAlgebra::norm2(r_local);
  });
    
  first_thread.join();
  second_thread.join();
}