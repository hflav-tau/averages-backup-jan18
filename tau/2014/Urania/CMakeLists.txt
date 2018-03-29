################################################################################
# Package: MCLimit
################################################################################
project(MCLimit)
cmake_minimum_required(VERSION 3.7)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations)

find_package(ROOT COMPONENTS Core RIO MathCore Minuit Hist Graf MathCore Minuit Hist Graf)

include(${ROOT_USE_FILE})
REFLEX_GENERATE_DICTIONARY(MCLimit dict/MCLimitDict.h SELECTION dict/MCLimit.xml)

include_directories(${CMAKE_CURRENT_LIST_DIR})

add_library(MCLimit SHARED
            src/junkInterface.cpp
            MCLimit.cxx)

target_link_libraries(MCLimit ${ROOT_LIBRARIES})

include(GNUInstallDirs)

install(TARGETS MCLimit DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY python DESTINATION ${CMAKE_INSTALL_PREFIX})
install(CODE "execute_process( \
         COMMAND ${CMAKE_COMMAND} -E create_symlink \
         libMCLimit.so \
         libMCLimitDict.so \
         )"
       )

install(FILES ${CMAKE_BINARY_DIR}/MCLimitDict.rootmap DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${CMAKE_BINARY_DIR}/MCLimit_rdict.pcm DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${CMAKE_BINARY_DIR}/libMCLimitDict.so DESTINATION ${CMAKE_INSTALL_LIBDIR})
