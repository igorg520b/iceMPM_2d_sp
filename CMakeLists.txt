cmake_minimum_required(VERSION 3.22)

project(WaveIce_MPM LANGUAGES CXX C CUDA)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored -use_fast_maths")
set(CMAKE_INCLUDE_CURRENT_DIR ON)

option(BUILD_CLI_VERSION "Build CLI version" OFF)
option(BUILD_GUI_VERSION "Build Qt version" ON)

if(BUILD_GUI_VERSION)
    set(CMAKE_AUTOUIC ON)
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTORCC ON)
endif()

#set(CUDA_TOOLKIT_ROOT_DIR "/usr/local/cuda")
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})




find_package(Eigen3 REQUIRED NO_MODULE HINTS "C:\\libraries\\eigen-build")
if(WIN32)
    add_compile_definitions(H5_BUILT_AS_DYNAMIC_LIB)
else()
    find_package(HDF5 REQUIRED COMPONENTS CXX)
    find_package(spdlog REQUIRED)
endif()



if(BUILD_CLI_VERSION)
# CLI version
add_executable(cm2
    cli/main.cpp

    simulation/parameters_sim.h
    simulation/parameters_sim.cpp
    simulation/model.cpp
    simulation/model.h
    simulation/point.h
    simulation/point.cpp

    simulation/givens.cuh
    simulation/helper_math.cuh
    simulation/gpu_implementation5.h
    simulation/gpu_implementation5.cpp

    simulation/snapshotmanager.h
    simulation/snapshotmanager.cpp
    simulation/gpu_partition.h
    simulation/gpu_partition.cpp
    simulation/kernels.cuh

    simulation/hssoa/host_side_soa.h
    simulation/hssoa/host_side_soa.cpp
    simulation/hssoa/proxypoint2d.h
    simulation/hssoa/proxypoint2d.cpp
)
set_source_files_properties(simulation/gpu_partition.cpp PROPERTIES LANGUAGE CUDA)

set_target_properties(cm2 PROPERTIES CUDA_ARCHITECTURES "80;89")
set_target_properties(cm2 PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_RESOLVE_DEVICE_SYMBOLS ON)

# spdlog
if (WIN32)
    target_link_libraries(cm2 PRIVATE spdlog)
    target_include_directories(cm2 PRIVATE "C:\\libraries\\spdlog\\include")
    target_link_directories(cm2 PRIVATE "C:\\libraries\\spdlog-build\\Release")
else()
    target_link_libraries(cm2 PRIVATE spdlog::spdlog)
endif()

target_include_directories(cm2 PRIVATE simulation simulation/hssoa)
target_link_libraries(cm2 PRIVATE Eigen3::Eigen)

# HDF5
if(WIN32)
    target_link_directories(cm2 PRIVATE "C:\\libraries\\hdf5-build\\bin\\Release")
    target_link_libraries(cm2 PRIVATE hdf5_cpp hdf5 )
else()
    target_link_libraries(cm2 PRIVATE ${HDF5_LIBRARIES})
    target_include_directories(cm2 PRIVATE ${HDF5_INCLUDE_DIRS})
endif()

endif()


if(BUILD_GUI_VERSION)
# Qt GUI executable

include(FetchContent)
FetchContent_Declare(
    PBrowser
    GIT_REPOSITORY https://github.com/igorg520b/PBrowser.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(PBrowser)


    find_package(VTK 9.1 HINTS "C:\\libraries\\vtk-build\\lib\\cmake\\vtk-9.1" REQUIRED COMPONENTS
      CommonColor
      CommonCore
      FiltersSources
      FiltersGeometry
      InteractionStyle
      RenderingContextOpenGL2
      RenderingCore
      RenderingFreeType
      RenderingGL2PSOpenGL2
      RenderingOpenGL2
      IOGeometry
      RenderingAnnotation
      RenderingLabel
      IOImage
      ViewsQt
    )
    find_package(Qt5 COMPONENTS Widgets REQUIRED)
    add_compile_definitions(QT_NO_KEYWORDS)

    add_executable(gm2
        gui/main.cpp

        gui/windows/mainwindow.cpp
        gui/windows/mainwindow.h
        gui/windows/mainwindow.ui
        icons.qrc

        gui/vtk/vtk_representation.cpp
        gui/vtk/vtk_representation.h

        gui/backgroundworker.cpp
        gui/backgroundworker.h
        gui/parameters_wrapper.h


        # simulation
        simulation/parameters_sim.h
        simulation/parameters_sim.cpp
        simulation/model.cpp
        simulation/model.h
        simulation/point.h
        simulation/point.cpp

        simulation/givens.cuh
        simulation/helper_math.cuh
        simulation/gpu_implementation5.h
        simulation/gpu_implementation5.cpp

        simulation/snapshotmanager.h
        simulation/snapshotmanager.cpp

        simulation/hssoa/host_side_soa.h
        simulation/hssoa/host_side_soa.cpp
        simulation/hssoa/proxypoint2d.h
        simulation/hssoa/proxypoint2d.cpp

        simulation/gpu_partition.h
        simulation/gpu_partition.cpp
        simulation/kernels.cuh

        simulation/poisson_disk_sampling.h

    )
    set_source_files_properties(simulation/gpu_partition.cpp PROPERTIES LANGUAGE CUDA)
    
    message(${PBrowser_SOURCE_DIR})
    target_include_directories(gm2 PRIVATE ${PBrowser_SOURCE_DIR} ${PBrowser_SOURCE_DIR}/propertybrowser)
    target_include_directories(gm2 PRIVATE ${PBrowser_SOURCE_DIR}/objectpropertybrowser)
    target_link_libraries(gm2 PRIVATE PBrowser)

    set_target_properties(gm2 PROPERTIES CUDA_ARCHITECTURES "80;89")
    set_target_properties(gm2 PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_RESOLVE_DEVICE_SYMBOLS ON)

    target_include_directories(gm2 PRIVATE gui gui/windows gui/vtk simulation simulation/hssoa External)
    target_link_libraries(gm2 PRIVATE Qt5::Widgets Qt5::Core)

    # Eigen
    find_package(Eigen3 REQUIRED NO_MODULE HINTS "C:\\libraries\\eigen-build")
    target_link_libraries(gm2 PRIVATE Eigen3::Eigen)

    # VTK
    target_link_libraries(gm2 PRIVATE ${VTK_LIBRARIES})



    # HDF5
    if(WIN32)
        target_link_directories(gm2 PRIVATE "C:\\libraries\\hdf5-build\\bin\\Release")
        target_link_libraries(gm2 PRIVATE hdf5_cpp hdf5 )
    else()
        target_link_libraries(gm2 PRIVATE ${HDF5_LIBRARIES})
        target_include_directories(gm2 PRIVATE ${HDF5_INCLUDE_DIRS})
    endif()

    # spdlog
    if (WIN32)
        target_link_libraries(gm2 PRIVATE spdlog)
        target_include_directories(gm2 PRIVATE "C:\\libraries\\spdlog\\include")
        target_link_directories(gm2 PRIVATE "C:\\libraries\\spdlog-build\\Release")
    else()
        target_link_libraries(gm2 PRIVATE spdlog::spdlog)
    endif()


    # OpenMP
#[[    if(WIN32)
        target_compile_options(iceMPM_GUI PRIVATE /openmp)
    else()
        target_compile_options(iceMPM_GUI PRIVATE -fopenmp -Wno-stringop-overflow -Wno-invalid-offsetof)
        target_link_libraries(iceMPM_GUI PRIVATE gomp)
    endif()
]]
vtk_module_autoinit(
  TARGETS gm2
  MODULES ${VTK_LIBRARIES}
)

endif()

