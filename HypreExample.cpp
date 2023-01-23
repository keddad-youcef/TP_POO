#include <iostream>
#include "HypreExample.h"
#include "AlienMock.h"
#include "LocalLinearAlgebra.h"
#include "Timer.h"


//TD2 & TD3//


ResidualNorms HypreExample::run()
{

  // Start counting for total calcul time...
  Timer timer_global{};
  std::chrono::time_point<std::chrono::system_clock> start_global, end_global;
  start_global = timer_global.start();
  //


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

  // dimensionnement et initialisation à 0 de la matrice
  LocalLinearAlgebra::Matrix A_local(size, size, 0);
  

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
      A_local.add_value(irow, irow,2.);
      if (irow - 1 >= 0)
        builder(irow, irow - 1) = -1.;
        A_local.add_value(irow, irow - 1,-1.);
      if (irow + 1 < gsize)
        builder(irow, irow + 1) = -1.;
        A_local.add_value(irow, irow + 1,-1.);
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



  // Start counting for locallinearalgebra calcul time...
  Timer timer_local{};
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = timer_local.start();
  //


  // LocalLinearAlgebra implementation

  LocalLinearAlgebra::Vector tmp_local(size);
  
  LocalLinearAlgebra::Vector r_local;

  LocalLinearAlgebra::Vector b_local(size);

  // dimensionnement et initialisation à 1 du vecteur
  LocalLinearAlgebra::Vector x_local(size, 1);
  
 

  LocalLinearAlgebra::mult(A_local, x_local, tmp_local); 
  
  LocalLinearAlgebra::mult(A_local, x_local, b_local);

  
  //
  LocalLinearAlgebra::copy(tmp_local, r_local);

  //
  LocalLinearAlgebra::axpy(-1., b_local, r_local);



  auto norm_local = LocalLinearAlgebra::norm2(r_local);


  //
  end = timer_local.stop();          // Stop counting for locallinearalgebra calcul time...
  std::cout << "LocalLinearAlgebra ";
  timer_local.print(end, start);     // Print the time spent. 

  //

  end_global = timer_global.stop();   // Stop counting for total calcul time...
  std::cout << "Total ";
  timer_global.print(end_global, start_global);   // Print the time spent.


  std::cout << "\n\n\n";
  
  ResidualNorms ResNormObj{norm,norm_local};

  return ResNormObj;
};







