if(NOT DEFINED ENV{CXXFLAGS})
    if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unknown-pragmas -Wno-sign-compare -Woverloaded-virtual -Wwrite-strings -Wno-unused -fdiagnostics-color=always")
        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
        set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g3 -march=native -DNDEBUG")
        set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -DDEBUG")
    endif()
endif()
