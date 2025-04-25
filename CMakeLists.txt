cmake_minimum_required(VERSION 3.20.2)

project(WaveIce_MPM LANGUAGES CXX C CUDA)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored -use_fast_maths")
set(CMAKE_INCLUDE_CURRENT_DIR ON)

option(BUILD_CLI_VERSION "Build CLI version" ON)
option(BUILD_GUI_VERSION "Build Qt version" OFF)
option(BUILD_POSTPROCESSOR "Build Qt version" OFF)

if(BUILD_GUI_VERSION OR BUILD_POSTPROCESSOR)
    set(CMAKE_AUTOUIC ON)
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTORCC ON)
endif()

include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(spdlog REQUIRED)
find_package(fmt REQUIRED)



if(BUILD_CLI_VERSION)

# CLI version
add_executable(cm2
    cli/main.cpp

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


    simulation/hssoa/host_side_soa.h
    simulation/hssoa/host_side_soa.cpp
    simulation/hssoa/proxypoint2d.h
    simulation/hssoa/proxypoint2d.cpp

    simulation/gpu_partition.h
    simulation/gpu_partition.cpp
    simulation/kernels.cuh

    simulation/data_manager/snapshotmanager.h
    simulation/data_manager/snapshotmanager.cpp
    simulation/data_manager/windandcurrentinterpolator.cpp
    simulation/data_manager/windandcurrentinterpolator.h

    simulation/data_manager/poisson_disk_sampling.h
)
set_source_files_properties(simulation/gpu_partition.cpp PROPERTIES LANGUAGE CUDA)

set_target_properties(cm2 PROPERTIES CUDA_ARCHITECTURES "80;89")
set_target_properties(cm2 PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_RESOLVE_DEVICE_SYMBOLS ON)

target_link_libraries(cm2 PRIVATE Eigen3::Eigen spdlog::spdlog)
target_link_libraries(cm2 PRIVATE fmt::fmt)

target_include_directories(cm2 PRIVATE simulation simulation/hssoa simulation/data_manager External External/rapidjson)


if(HDF5_INCLUDE_DIRS)
    target_link_libraries(cm2 PRIVATE ${HDF5_LIBRARIES})
else()
#    target_link_libraries(cm2 PRIVATE hdf5-shared hdf5_cpp-shared)
    target_link_libraries(cm2 PRIVATE hdf5_cpp-shared)
endif()

target_include_directories(cm2 PRIVATE ${HDF5_INCLUDE_DIRS})

endif()










if(BUILD_POSTPROCESSOR)
    find_package(VTK 9.4 HINTS "~/Projects-libraries/build-vtk/lib/cmake/vtk-9.4" REQUIRED COMPONENTS
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

    add_executable(pp
    	postprocessor/pp_main.cpp
    	postprocessor/window/pp_mainwindow.cpp
    	postprocessor/window/pp_mainwindow.h
    	postprocessor/window/pp_mainwindow.ui
        icons.qrc
        
        postprocessor/vtk/vtk_visualization.cpp
        postprocessor/vtk/vtk_visualization.h

        # simulation
        simulation/parameters_sim.h
        simulation/parameters_sim.cpp

        postprocessor/framedata.h
        postprocessor/framedata.cpp
        postprocessor/generalgriddata.h
        postprocessor/generalgriddata.cpp

        gui/colormap.h
        gui/colormap.cpp
    )

    target_include_directories(pp PRIVATE postprocessor postprocessor/window postprocessor/vtk gui External)
    target_include_directories(pp PRIVATE simulation simulation/hssoa simulation/data_manager)
    target_link_libraries(pp PRIVATE Qt5::Widgets Qt5::Core)

    # Eigen
    find_package(Eigen3 REQUIRED NO_MODULE)
    target_link_libraries(pp PRIVATE Eigen3::Eigen)

    # VTK
    target_link_libraries(pp PRIVATE ${VTK_LIBRARIES})

    # HDF5
    target_link_libraries(pp PRIVATE ${HDF5_LIBRARIES})
    target_include_directories(pp PRIVATE ${HDF5_INCLUDE_DIRS})

    # spdlog
    target_link_libraries(pp PRIVATE spdlog::spdlog)

vtk_module_autoinit(
  TARGETS pp
  MODULES ${VTK_LIBRARIES}
)

target_compile_options(pp PRIVATE -fopenmp -Wno-stringop-overflow -Wno-invalid-offsetof)
target_link_libraries(pp PRIVATE gomp)

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


#    find_package(VTK 9.4 HINTS "~/Projects-libraries/build-vtk/lib/cmake/vtk-9.4" REQUIRED COMPONENTS
        find_package(VTK 9.1 REQUIRED COMPONENTS
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


        simulation/hssoa/host_side_soa.h
        simulation/hssoa/host_side_soa.cpp
        simulation/hssoa/proxypoint2d.h
        simulation/hssoa/proxypoint2d.cpp

        simulation/gpu_partition.h
        simulation/gpu_partition.cpp
        simulation/kernels.cuh
        
        simulation/data_manager/snapshotmanager.h
        simulation/data_manager/snapshotmanager.cpp
        simulation/data_manager/windandcurrentinterpolator.cpp
        simulation/data_manager/windandcurrentinterpolator.h
        simulation/data_manager/poisson_disk_sampling.h


        gui/vtk/SpecialSelector2D.cpp
        gui/vtk/SpecialSelector2D.h
        gui/colormap.h
        gui/colormap.cpp
    )
    set_source_files_properties(simulation/gpu_partition.cpp PROPERTIES LANGUAGE CUDA)
    
    message(${PBrowser_SOURCE_DIR})
    target_include_directories(gm2 PRIVATE ${PBrowser_SOURCE_DIR} ${PBrowser_SOURCE_DIR}/propertybrowser)
    target_include_directories(gm2 PRIVATE ${PBrowser_SOURCE_DIR}/objectpropertybrowser)
    target_link_libraries(gm2 PRIVATE PBrowser)

    set_target_properties(gm2 PROPERTIES CUDA_ARCHITECTURES "80;89")
    set_target_properties(gm2 PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_RESOLVE_DEVICE_SYMBOLS ON)

    target_include_directories(gm2 PRIVATE gui gui/windows gui/vtk simulation simulation/hssoa simulation/data_manager External External/rapidjson)
    target_link_libraries(gm2 PRIVATE Qt5::Widgets Qt5::Core)

    # Eigen
    find_package(Eigen3 REQUIRED NO_MODULE)
    target_link_libraries(gm2 PRIVATE Eigen3::Eigen)

    # VTK
    target_link_libraries(gm2 PRIVATE ${VTK_LIBRARIES})



    # HDF5
    target_link_libraries(gm2 PRIVATE ${HDF5_LIBRARIES})
    target_include_directories(gm2 PRIVATE ${HDF5_INCLUDE_DIRS})

    # spdlog
    target_link_libraries(gm2 PRIVATE spdlog::spdlog)


    # OpenMP
#[[
        target_compile_options(iceMPM_GUI PRIVATE -fopenmp -Wno-stringop-overflow -Wno-invalid-offsetof)
        target_link_libraries(iceMPM_GUI PRIVATE gomp)
]]
vtk_module_autoinit(
  TARGETS gm2
  MODULES ${VTK_LIBRARIES}
)

endif()

