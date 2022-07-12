# OpenDataDetector

[![](https://zenodo.org/badge/DOI/10.5281/zenodo.4674401.svg)](https://doi.org/10.5281/zenodo.4674401)

The `OpenDataDetector` (ODD) is attempted to provide a template (HL-)LHC style particle detector for algorithm research and development.

## Tracking Detector

The ODD Tracking system is an evolution of the detector used for the `Tracking Machine Learning Challenge` (part 1 and 2), and comprises a system of several components:
 * a central beam pipe
 * an innermost Pixel tracking system
 * a middle Short Strip system
 * an outermost Long Strip system
 * an enclosing solenoid 

 ## Build instructions

 The ODD library can be built using `CMake` with minimal dependencies (mainly required by DD4hep), dependencies are:
 * BOOST
 * DD4hep
 * ROOT
 * Geant4

 ### Building with CMake    

The following will build the ODD DD4hep detector:

```shell
cmake -S <path_to_source> -B <path_to_build_area>  -DDD4hep_DIR=<path_to_DD4hp> cmake -DGeant4_DIR=<path_to_Geant4> -DROOT_DIR=<path_to_ROOT> -DCMAKE_CXX_STANDARD=17
cmake --build <path_to_build_area>
 ```


