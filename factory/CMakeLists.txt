add_library(
    OpenDataDetector SHARED
    tracker/ODDCylinder_geo.cpp
    tracker/ODDModuleHelper.cpp
    tracker/ODDPixelBarrel_geo.cpp
    tracker/ODDPixelEndcap_geo.cpp
    tracker/ODDStripBarrel_geo.cpp
    tracker/ODDStripEndcap_geo.cpp)

    
target_link_libraries(OpenDataDetector PUBLIC DD4hep::DDCore DD4hep::DDRec)

# Not sure why this needs to be set, but dd4hep fails to compile otherwise
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})
set_target_properties(OpenDataDetector PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
dd4hep_set_version(OpenDataDetector MAJOR 1 MINOR 0 PATCH 0)
dd4hep_generate_rootmap(OpenDataDetector)

install(
  TARGETS OpenDataDetector
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
