cmake_minimum_required(VERSION 2.6)

# Source code here
add_sources(brODErs++
    DataShocking1T.cpp
    DataShockingNT.cpp
    Mesh.cpp
    RankineHugoniot1T.cpp
    RankineHugoniotNT.cpp
    SetupShocking1T.cpp
    SetupShockingNT.cpp
    Shocking1T.cpp
    ShockingTTv.cpp
)

install(FILES DataShocking1T.h DESTINATION bin/include)
install(FILES Mesh.h DESTINATION bin/include)
install(FILES Shocking1T.h DESTINATION bin/include)
install(FILES ShockingTTv.h DESTINATION bin/include)
install(FILES ShockRelations.h DESTINATION bin/include)
install(FILES RankineHugoniot1T.h DESTINATION bin/include)
install(FILES SetupShocking1T.h DESTINATION bin/include)
install(FILES SetupShockingNT.h DESTINATION bin/include)

