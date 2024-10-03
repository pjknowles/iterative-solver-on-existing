cmake_minimum_required(VERSION 3.28)
project(iterative_solver_on_existing Fortran CXX)
include(FetchContent)
FetchContent_Declare(iterative-solver
        GIT_REPOSITORY https://github.com/knowles-group/iterative-solver
        GIT_TAG one-call-interface)
FetchContent_MakeAvailable(iterative-solver)

enable_language(Fortran)
add_executable(iterative_solver_on_existing main.f90)
target_compile_options(iterative_solver_on_existing PRIVATE -fdefault-integer-8)
target_link_libraries(iterative_solver_on_existing PRIVATE iterative-solver)