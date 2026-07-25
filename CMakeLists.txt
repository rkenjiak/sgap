cmake_minimum_required(VERSION 3.20)
project(sgap LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# =========================================================================
#  sgap-core: biblioteca estatica (core + algorithms + io)
# =========================================================================
add_library(sgap-core STATIC
    core/sequence.cpp
    core/sequence_graph.cpp
    core/alignment_graph.cpp
    algorithms/alignment_walk.cpp
    algorithms/alignment_potentials.cpp
    algorithms/warm_start.cpp
    io/gfa.cpp
    io/fastq.cpp
    io/config.cpp
    io/gurobi_sol.cpp
)
target_include_directories(sgap-core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# =========================================================================
#  sgap-heuristic: executavel de heuristicas (sem dependencias externas)
# =========================================================================
add_executable(sgap-heuristic
    apps/main.cpp
    heuristics/greedy_degree.cpp
    heuristics/beam_search.cpp
    # Relaxacao Lagrangeana: subgradiente + Dijkstra + beam de reparo.
    # Nao usa Gurobi, entao vive aqui (nao no sgap-solver) e fica
    # disponivel sem BUILD_SOLVER=ON.
    solvers/lagrangian.cpp
)
target_link_libraries(sgap-heuristic PRIVATE sgap-core)

# =========================================================================
#  sgap-solver: executavel PLI + relaxacoes (requer Gurobi, opcional)
# =========================================================================
option(BUILD_SOLVER "Build exact solver — requires Gurobi" OFF)

if(BUILD_SOLVER)
    # -----------------------------------------------------------------
    # Deteccao do Gurobi
    # -----------------------------------------------------------------
    if(NOT DEFINED GUROBI_HOME)
        if(DEFINED ENV{GUROBI_HOME})
            set(GUROBI_HOME "$ENV{GUROBI_HOME}")
        else()
            file(GLOB _gurobi_candidates "/opt/gurobi*/linux64")
            if(_gurobi_candidates)
                list(SORT _gurobi_candidates ORDER DESCENDING)
                list(GET _gurobi_candidates 0 GUROBI_HOME)
            endif()
        endif()
    endif()

    if(NOT GUROBI_HOME)
        message(FATAL_ERROR
            "BUILD_SOLVER=ON mas Gurobi nao encontrado.\n"
            "Defina GUROBI_HOME, por exemplo:\n"
            "  cmake -B build -DBUILD_SOLVER=ON -DGUROBI_HOME=/opt/gurobi1301/linux64\n"
            "ou exporte a variavel de ambiente:\n"
            "  export GUROBI_HOME=/opt/gurobi1301/linux64")
    endif()

    message(STATUS "Usando Gurobi em: ${GUROBI_HOME}")

    find_path(GUROBI_INCLUDE_DIR
        NAMES gurobi_c++.h
        PATHS "${GUROBI_HOME}/include"
        NO_DEFAULT_PATH)
    if(NOT GUROBI_INCLUDE_DIR)
        message(FATAL_ERROR "gurobi_c++.h nao encontrado em ${GUROBI_HOME}/include")
    endif()

    file(GLOB _gurobi_c_libs "${GUROBI_HOME}/lib/libgurobi[0-9]*.so")
    if(NOT _gurobi_c_libs)
        message(FATAL_ERROR "libgurobiXXX.so nao encontrada em ${GUROBI_HOME}/lib")
    endif()
    list(SORT _gurobi_c_libs ORDER DESCENDING)
    list(GET _gurobi_c_libs 0 GUROBI_C_LIB)

    find_library(GUROBI_CXX_LIB
        NAMES gurobi_c++
        PATHS "${GUROBI_HOME}/lib"
        NO_DEFAULT_PATH)
    if(NOT GUROBI_CXX_LIB)
        message(FATAL_ERROR "libgurobi_c++.a nao encontrada em ${GUROBI_HOME}/lib")
    endif()

    message(STATUS "  Include : ${GUROBI_INCLUDE_DIR}")
    message(STATUS "  Lib C   : ${GUROBI_C_LIB}")
    message(STATUS "  Lib C++ : ${GUROBI_CXX_LIB}")

    # -----------------------------------------------------------------
    # Target
    # -----------------------------------------------------------------
    add_executable(sgap-solver
        apps/main-pli.cpp
        solvers/pli.cpp
        # solvers/lagrangian.cpp  -> movido para sgap-heuristic (sem Gurobi)
        # Descomentar conforme implementar:
        # solvers/column_gen.cpp
        # solvers/benders.cpp
    )
    target_include_directories(sgap-solver PRIVATE "${GUROBI_INCLUDE_DIR}")
    target_link_libraries(sgap-solver PRIVATE
        sgap-core
        "${GUROBI_CXX_LIB}"
        "${GUROBI_C_LIB}"
        pthread m dl
    )
    set_target_properties(sgap-solver PROPERTIES
        BUILD_RPATH  "${GUROBI_HOME}/lib"
        INSTALL_RPATH "${GUROBI_HOME}/lib"
    )

    message(STATUS "sgap-solver: ON")
else()
    message(STATUS "sgap-solver: OFF (use -DBUILD_SOLVER=ON para habilitar)")
endif()
