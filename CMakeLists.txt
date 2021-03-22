cmake_minimum_required(VERSION 3.9)
project(spssLangevinIntegration)

set(CMAKE_CXX_STANDARD 17)

include_directories(/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include
        /Users/nikita/Documents/mpich-install/include/ /usr/local/include/)
link_directories(/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/lib
        /Users/nikita/Documents/mpich-install/lib/ /usr/local/lib/)

set(GCC_COVERAGE_COMPILE_FLAGS "-march=native -O3 -flto")
set(GCC_COVERAGE_LINK_FLAGS "")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
set(CMAKE_OSX_SYSROOT /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk)

add_executable(spssLangevinIntegration
        main.cpp
        Definitions.hpp
        Observers/BinaryObserverSecondOrder.cpp Observers/BinaryObserverSecondOrder.hpp
        BoundaryConfiguration/PeriodicBoundaryConditions.cpp BoundaryConfiguration/PeriodicBoundaryConditions.hpp
        Steppers/RungeKutta4Stepper.cpp Steppers/RungeKutta4Stepper.hpp
        DynamicalSystems/ParticleSystemSecondOrder.cpp DynamicalSystems/ParticleSystemSecondOrder.hpp
        ControlEngines/SimulationEngineSecondOrder.cpp ControlEngines/SimulationEngineSecondOrder.hpp
        Parallelization/Thread.cpp Parallelization/Thread.hpp
        Parallelization/Parallelization.hpp Parallelization/Parallelization.cpp
        Parallelization/ThreadSharedMemory.cpp Parallelization/ThreadSharedMemory.hpp
        Observers/BinaryObserverSecondOrderPtr.cpp Observers/BinaryObserverSecondOrderPtr.hpp
        Steppers/RungeKutta4StepperPtr.cpp Steppers/RungeKutta4StepperPtr.hpp
        DynamicalSystems/ParticleSystemSecondOrderPtr.cpp DynamicalSystems/ParticleSystemSecondOrderPtr.hpp
        ControlEngines/SimulationEngineSecondOrderPtr.cpp ControlEngines/SimulationEngineSecondOrderPtr.hpp
        Steppers/StochasticEulerStepper.hpp
        Steppers/StochasticEulerStepperPtr.hpp
        ControlEngines/SimulationEngineFirstOrder.cpp ControlEngines/SimulationEngineFirstOrder.hpp
        DynamicalSystems/ParticleSystemFirstOrder.cpp DynamicalSystems/ParticleSystemFirstOrder.hpp
        Observers/BinaryObserverFirstOrder.cpp Observers/BinaryObserverFirstOrder.hpp
        ControlEngines/SimulationEngineFirstOrderPtr.cpp ControlEngines/SimulationEngineFirstOrderPtr.hpp
        DynamicalSystems/ParticleSystemFirstOrderPtr.cpp DynamicalSystems/ParticleSystemFirstOrderPtr.hpp
        Observers/BinaryObserverFirstOrderPtr.cpp Observers/BinaryObserverFirstOrderPtr.hpp
        Observers/AveragingObserver.cpp Observers/AveragingObserver.hpp
        Steppers/StochasticRungeKuttaStepperPtr.hpp
        Observers/AveragingObserverPtr.cpp Observers/AveragingObserverPtr.hpp)