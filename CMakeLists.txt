#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(G4PlasticDetector)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
# Link against OpenCV libraries
#
add_executable(G4PlasticDetector G4PlasticDetector.cc ${sources} ${headers})
target_link_libraries(G4PlasticDetector ${Geant4_LIBRARIES})

set(G4PLASTICDETECTOR_SCRIPTS
  macros/vis.mac
  macros/rotate.mac
  macros/runElectronUniformLinISO.mac
  macros/runElectronUniformLogISO.mac
  macros/runElectronPowerISO.mac
  macros/runGammaPowerElectronGaussianISO.mac
  macros/runGammaUniformLinISO.mac
  macros/runGammaUniformLogISO.mac
  macros/runGammaPowerISO.mac
  macros/runCl36ISO.mac
  macros/runCo60ISO.mac
  macros/runCs137ISO.mac
  macros/runSr90Y90ISO.mac
)

foreach(_script ${G4PLASTICDETECTOR_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS G4PlasticDetector DESTINATION bin )
