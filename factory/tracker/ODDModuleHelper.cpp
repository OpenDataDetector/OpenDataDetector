// Open Data Dector project
//
// (c) 2021 CERN for the benefit of the ODD project
//
// Mozilla Public License Version 2.0

#include "ODDModuleHelper.hpp"

#include "DDRec/DetectorData.h"
#include "ODDHelper.hpp"

using namespace dd4hep;

std::pair<Assembly, DetElement> ODDModuleHelper::assembleTrapezoidalModule(
    Detector &oddd, SensitiveDetector &sens, const xml_comp_t &x_module) {
  // The Module envelope volume
  Assembly moduleAssembly("module");
  // Visualization
  moduleAssembly.setVisAttributes(oddd, x_module.visStr());

  // The module detector element
  DetElement moduleElement("ModuleElementTemplate", 0);

  // Place the components inside the module
  unsigned int compNum = 0;
  unsigned int sensorNum = 0;

  for (xml_coll_t comp(x_module, _U(module_component)); comp;
       ++comp, ++compNum) {
    xml_comp_t x_comp = comp;

    // create the component volume
    std::string compName =
        _toString((int)compNum, "component%d") + x_comp.materialStr();

    Trapezoid trapShape(x_comp.x1(), x_comp.x2(), 0.5 * x_comp.thickness(),
                        0.5 * x_comp.thickness(), x_comp.length());

    Volume componentVolume(compName, trapShape,
                           oddd.material(x_comp.materialStr()));
    componentVolume.setVisAttributes(oddd, x_comp.visStr());

    // Place carbon foam structure
    if (x_comp.hasChild(_U(subtraction)) and x_comp.hasChild(_U(tube))) {
      xml_comp_t x_sub = x_comp.child(_U(subtraction));
      xml_comp_t x_tubs = x_sub.child(_U(tubs));
      xml_comp_t x_pipe = x_comp.child(_U(tube));
      double length = x_comp.length();

      // Create the subtraction
      Tube tubeCutoutSeg1(x_tubs.rmin(), x_tubs.rmax(), length + x_tubs.dz(),
                          -0.1, 2.1 * M_PI);
      Tube tubeCutoutSeg2(x_tubs.rmin(), x_tubs.rmax(), length + x_tubs.dz(), 0,
                          2 * M_PI);
      UnionSolid foamCutout(tubeCutoutSeg1, tubeCutoutSeg2);

      componentVolume = Volume(
          compName,
          SubtractionSolid(
              trapShape, foamCutout,
              Position(x_sub.x_offset(), x_sub.y_offset(), x_sub.z_offset())),
          oddd.material(x_comp.materialStr()));
      componentVolume.setVisAttributes(oddd, x_comp.visStr());

      // Create the cooling pipe
      Tube coolingPipe(x_pipe.rmin(), x_pipe.rmax(),
                       0.5 * length + x_pipe.dz());
      Volume pipeVolume("CoolingPipe", coolingPipe,
                        oddd.material(x_pipe.materialStr()));
      pipeVolume.setVisAttributes(oddd, x_pipe.visStr());

      // Place the pipe in the module
      moduleAssembly.placeVolume(
          pipeVolume,
          Position(x_pipe.x_offset(), x_pipe.y_offset(), x_pipe.z_offset()));
    }

    // Place the component
    double stereoAlpha = x_comp.alpha();
    PlacedVolume placedComponent = moduleAssembly.placeVolume(
        componentVolume,
        Transform3D(
            RotationY(stereoAlpha),
            Position(x_comp.x_offset(), x_comp.y_offset(), x_comp.z_offset())));

    // Deal with the sensitive sensor
    if (x_comp.isSensitive()) {
      sens.setType("tracker");
      componentVolume.setSensitiveDetector(sens);
      placedComponent.addPhysVolID("sensor", sensorNum++);
      // Create the sensor element and place it
      std::string sensorName = _toString((int)sensorNum, "sensor%d");
      DetElement sensorElement(moduleElement, sensorName, sensorNum);
      sensorElement.setPlacement(placedComponent);

      // Add the sensor extension
      auto &params = ODDHelper::ensureExtension<dd4hep::rec::VariantParameters>(
          sensorElement);
      params.set<std::string>("axis_definitions", "XZY");
    }
  }

  // return the module assembly
  return std::pair<Assembly, DetElement>(moduleAssembly, moduleElement);
}

std::pair<Assembly, DetElement> ODDModuleHelper::assembleRectangularModule(
    Detector &oddd, SensitiveDetector &sens, const xml_comp_t &x_module,
    double &ylength) {
  // The Module envelope volume
  Assembly moduleAssembly("module");
  // Visualization
  moduleAssembly.setVisAttributes(oddd, x_module.visStr());

  // The module detector element
  DetElement moduleElement("ModuleElementTemplate", 0);

  // Place the components inside the module
  unsigned int compNum = 0;
  unsigned int sensorNum = 0;

  for (xml_coll_t comp(x_module, _U(module_component)); comp;
       ++comp, ++compNum) {
    xml_comp_t x_comp = comp;

    // Component volume
    std::string componentName = _toString((int)compNum, "component%d");
    Box boxShape(0.5 * x_comp.dx(), 0.5 * x_comp.dy(), 0.5 * x_comp.dz());
    // Standard component volume without cutout
    Volume componentVolume(componentName, boxShape,
                           oddd.material(x_comp.materialStr()));

    // Place carbon foam structure
    if (x_comp.hasChild(_U(subtraction)) and x_comp.hasChild(_U(tube))) {
      xml_comp_t x_sub = x_comp.child(_U(subtraction));
      xml_comp_t x_tubs = x_sub.child(_U(tubs));
      xml_comp_t x_pipe = x_comp.child(_U(tube));
      double length = x_comp.dy();

      // Create the subtraction
      Tube tubeCutoutSeg1(x_tubs.rmin(), x_tubs.rmax(), length + x_tubs.dz(),
                          -0.1, 2.1 * M_PI);
      Tube tubeCutoutSeg2(x_tubs.rmin(), x_tubs.rmax(), length + x_tubs.dz(), 0,
                          2 * M_PI);
      UnionSolid foamCutout(tubeCutoutSeg1, tubeCutoutSeg2);

      componentVolume =
          Volume(componentName,
                 SubtractionSolid(
                     boxShape, foamCutout,
                     Transform3D(RotationX(0.5 * M_PI),
                                 Position(x_sub.x_offset(), x_sub.y_offset(),
                                          x_sub.z_offset()))),
                 oddd.material(x_comp.materialStr()));
      componentVolume.setVisAttributes(oddd, x_comp.visStr());

      // Create the cooling pipe
      Tube coolingPipe(x_pipe.rmin(), x_pipe.rmax(),
                       0.5 * length + x_pipe.dz());
      Volume pipeVolume("CoolingPipe", coolingPipe,
                        oddd.material(x_pipe.materialStr()));
      pipeVolume.setVisAttributes(oddd, x_pipe.visStr());

      // Place the pipe in the module
      moduleAssembly.placeVolume(
          pipeVolume, Transform3D(RotationX(0.5 * M_PI),
                                  Position(x_pipe.x_offset(), x_pipe.y_offset(),
                                           x_pipe.z_offset())));
    }

    componentVolume.setVisAttributes(oddd, x_comp.visStr());

    // Calculate the module dimension
    double cylength =
        2. * abs(std::copysign(0.5 * x_comp.dy(), x_comp.y_offset()) +
                 x_comp.y_offset());
    ylength = cylength > ylength ? cylength : ylength;

    // Visualization
    componentVolume.setVisAttributes(oddd, x_comp.visStr());
    // Place Module Box Volumes in layer
    double stereoAlpha = x_comp.alpha();
    PlacedVolume placedComponent = moduleAssembly.placeVolume(
        componentVolume,
        Transform3D(
            RotationZ(stereoAlpha),
            Position(x_comp.x_offset(), x_comp.y_offset(), x_comp.z_offset())));

    // Deal with the sensitive sensor
    if (x_comp.isSensitive()) {
      sens.setType("tracker");
      componentVolume.setSensitiveDetector(sens);
      placedComponent.addPhysVolID("sensor", sensorNum++);

      // Create the sensor element and place it
      std::string sensorName = _toString((int)sensorNum, "sensor%d");
      DetElement sensorElement(moduleElement, sensorName, sensorNum);
      sensorElement.setPlacement(placedComponent);

      // Add the sensor extension
      auto &params = ODDHelper::ensureExtension<dd4hep::rec::VariantParameters>(
          sensorElement);
      params.set<std::string>("axis_definitions", "XYZ");
    }
  }

  // return the module assembly
  return std::pair<Assembly, DetElement>(moduleAssembly, moduleElement);
}
