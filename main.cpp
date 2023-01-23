#include <iostream>
#include "HypreExample.h"
#include "PETScExample.h"
#include "MpiMock.h"
#include "GenericExample.h"

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

    //TD3

    HypreExample HypreExample{};
    PETScExample PETScExample{};
    
    auto HypRet = HypreExample.run();
    auto PetRet = PETScExample.run();



    // TD4
    GenericExample generic_example{};

    generic_example.info<HypreAPI>();
    auto GenericHyp = generic_example.run(SolverType::Hypre);

    generic_example.info<PETScAPI>();
    auto GenericPet = generic_example.run(SolverType::PETSc);
    


    // TD5
    generic_example.run_parallel_thread();
  

  MPI_Finalize();
  return 0;
  // return HypRet,PetRet;  
}

