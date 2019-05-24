cmake_minimum_required(VERSION 2.6)

# Source code here
add_sources(brODErs++
    BasicOutput.cpp
    main.cpp
    Solver.cpp
    SetupProblem.cpp
    SetupProperties.cpp
)

install(FILES BasicOutput.h DESTINATION bin/include)
install(FILES Data.h DESTINATION bin/include)
install(FILES dfdt.h DESTINATION bin/include)
install(FILES jacobian.h DESTINATION bin/include)
install(FILES Problem.h DESTINATION bin/include)
install(FILES SetupProblem.h DESTINATION bin/include)
install(FILES SetupProperties.h DESTINATION bin/include)
install(FILES Solver.h DESTINATION bin/include)
install(FILES Solution.h DESTINATION bin/include)

