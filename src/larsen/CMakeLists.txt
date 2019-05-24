cmake_minimum_required(VERSION 2.6)

# Source code here
add_sources(brODErs++
    DataLarsen1T.cpp
    DataLarsenTTv.cpp
    Larsen1T.cpp
    LarsenTTv.cpp
    SetupLarsen1T.cpp
    SetupLarsenTTv.cpp
)

install(FILES DataLarsen1T.h DESTINATION bin/include)
install(FILES DataLarsenTTv.h DESTINATION bin/include)
install(FILES Larsen1T.h DESTINATION bin/include)
install(FILES LarsenTTv.h DESTINATION bin/include)
install(FILES SetupLarsen1T.h DESTINATION bin/include)
install(FILES SetupLarsenTTv.h DESTINATION bin/include)

