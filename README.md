# OpenDataDetector

[![](https://zenodo.org/badge/DOI/10.5281/zenodo.4674401.svg)](https://doi.org/10.5281/zenodo.4674401)

The `OpenDataDetector` (ODD) is attempted to provide a template (HL-)LHC style particle detector for algorithm research and development.

## Sub-detector layout

The detector description is organized into staged XML fragments in `xml/`:

- `OpenDataDetector.xml`: full detector steering file (backward-compatible full build)
- `OpenDataDetectorDefs.xml`: global definitions bundle
  - material includes (`OpenDataDetectorElements.xml`, `OpenDataDetectorMaterials.xml`)
  - world and field constants
  - envelopes, identifiers, visualization
  - ACTS support include
- `OpenDataDetectorActsSupport.xml`: ACTS-specific support (material binning constants)
- `OpenDataDetectorTracker.xml`: tracker stage
  - includes `xml/detectors/BeamPipe.xml`
  - includes `xml/detectors/TrackerPixels.xml`
  - includes `xml/detectors/TrackerShortStrips.xml`
  - includes `xml/detectors/TrackerLongStrips.xml`
  - includes `xml/detectors/Solenoid.xml`
- `OpenDataDetectorCalorimeter.xml`: calorimeter stage
  - includes `xml/detectors/CalorimeterECal.xml`
  - includes `xml/detectors/CalorimeterHCal.xml`
- `OpenDataDetectorMuonSystem.xml`: muon stage
  - includes `xml/detectors/MuonSystem.xml`

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

### Displaying with DD4hep

You can use DD4hep geometry tools to load either the full detector or staged detector compositions.


### Full detector (backward compatible)

```sh
geoDisplay xml/OpenDataDetector.xml -load
```

### Staged detector loading

Tracker and solenoid only:

```sh
geoDisplay -input xml/OpenDataDetectorDefs.xml -input xml/OpenDataDetectorTracker.xml -load
```

Tracker, solenoid, and calorimeter:

```sh
geoDisplay -input xml/OpenDataDetectorDefs.xml -input xml/OpenDataDetectorTracker.xml -input xml/OpenDataDetectorCalorimeter.xml -load
```

Full staged chain (global + tracker + calo + muon):

```sh
geoDisplay -input xml/OpenDataDetectorDefs.xml -input xml/OpenDataDetectorTracker.xml -input xml/OpenDataDetectorCalorimeter.xml -input xml/OpenDataDetectorMuonSystem.xml -load
```

You can also display with `geoPluginRun`:

```sh
geoPluginRun -input xml/OpenDataDetector.xml  -interactive -plugin DD4hep_GeometryDisplay -level 8
```
