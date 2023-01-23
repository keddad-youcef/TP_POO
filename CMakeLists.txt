cmake_minimum_required(VERSION 3.10)

# set the project name
project(LinearAlgebraAPI)

# enable modern C++
set(CMAKE_CXX_STANDARD 17)




# add the executables
add_executable(main main.cpp HypreExample.cpp HypreExample.h PETScExample.cpp PETScExample.h LocalLinearAlgebra.cpp LocalLinearAlgebra.h GenericExample.cpp GenericExample.h Timer.cpp Timer.h)
add_executable(Test.exe Test.cpp HypreExample.cpp HypreExample.h PETScExample.cpp PETScExample.h LocalLinearAlgebra.cpp LocalLinearAlgebra.h GenericExample.cpp GenericExample.h Timer.cpp Timer.h)



#GTest Dependencies
find_package(GTest REQUIRED)
target_link_libraries(Test.exe GTest::GTest GTest::Main)


#Alien Dependencies (not used in my case only "AlienMock.h")
find_package(Alien REQUIRED)
target_link_libraries(main Alien::alien_semantic_ref Alien::alien_semantic_move Alien::hypre_wrapper Alien::petsc_wrapper)
target_link_libraries(Test.exe Alien::alien_semantic_ref Alien::alien_semantic_move Alien::hypre_wrapper Alien::petsc_wrapper)



#Testing Part
enable_testing()
add_test(Test Test.exe)
gtest_discover_tests(Test.exe)
