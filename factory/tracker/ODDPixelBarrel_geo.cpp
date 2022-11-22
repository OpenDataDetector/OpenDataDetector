// Open Data Dector project
//
// (c) 2021 CERN for the benefit of the ODD project
//
// Mozilla Public License Version 2.0

#include <vector>

#include "DD4hep/DetFactoryHelper.h"
#include "DDRec/DetectorData.h"
#include "ODDHelper.hpp"
#include "ODDModuleHelper.hpp"
#include "ODDServiceHelper.hpp"
#include "XML/Utilities.h"

using namespace std;
using namespace dd4hep;

/// Helper method to complete the stave structure with
/// cooling, support material, etc.
///
/// @param oddd the detector instance
/// @param x_stave the stave xml component
/// @param staveAssembly is the assembly to be completed
/// @param staveHlength the half length of the stave structure
/// @param ylength the length in y-direction
static void completeStaveStructure(Detector &oddd, xml_comp_t &x_stave,
                                   Assembly &staveAssembly, double staveHlength,
                                   double ylength) {
  unsigned int nModules = x_stave.nmodules();

  // Place carbon foam structure
  if (x_stave.hasChild(_U(subtraction)) and x_stave.hasChild(_U(tube))) {
    xml_comp_t x_sub = x_stave.child(_U(subtraction));
    xml_comp_t x_trd = x_sub.child(_U(trd));
    xml_comp_t x_tubs = x_sub.child(_U(tubs));
    xml_comp_t x_pipe = x_stave.child(_U(tube));

    // Create the two shapes first
    Trapezoid foamShape(x_trd.x1(), x_trd.x2(), staveHlength + x_trd.dz(),
                        staveHlength + x_trd.dz(), x_trd.thickness());

    // Create the subtraction
    Tube tubeCutoutSeg1(x_tubs.rmin(), x_tubs.rmax(),
                        staveHlength + x_tubs.dz(), -0.1, 2.1 * M_PI);
    Tube tubeCutoutSeg2(x_tubs.rmin(), x_tubs.rmax(),
                        staveHlength + x_tubs.dz(), 0, 2 * M_PI);
    UnionSolid foamCutout(tubeCutoutSeg1, tubeCutoutSeg2);

    // Create the subtraction
    Volume foamVolume("CarbonFoam",
                      SubtractionSolid(foamShape, foamCutout,
                                       Transform3D(RotationX(0.5 * M_PI))),
                      oddd.material(x_sub.materialStr()));
    foamVolume.setVisAttributes(oddd, x_sub.visStr());
    // Place the volume in the stave
    staveAssembly.placeVolume(
        foamVolume,
        Position(x_sub.x_offset(), x_sub.y_offset(), x_sub.z_offset()));

    Tube coolingPipe(x_pipe.rmin(), x_pipe.rmax(), staveHlength + x_pipe.dz());
    // Create the coolingpipe
    Volume pipeVolume("CoolingPipe", coolingPipe,
                      oddd.material(x_pipe.materialStr()));
    pipeVolume.setVisAttributes(oddd, x_pipe.visStr());

    // Place the pipe in the stave
    staveAssembly.placeVolume(
        pipeVolume, Transform3D(RotationX(0.5 * M_PI),
                                Position(x_pipe.x_offset(), x_pipe.y_offset(),
                                         x_pipe.z_offset())));

    xml_comp_t x_cable = x_stave.child(_U(eltube));
    // Place the support cables for the modules
    for (unsigned int modCable = 0; modCable < 0.5 * nModules; ++modCable) {
      double cableLength = staveHlength - modCable * ylength;

      Tube cable(x_cable.rmin(), x_cable.rmax(), 0.5 * cableLength);
      // Create the cable volume
      Volume cableVolume("Cable", cable,
                         oddd.material(x_cable.materialStr()));
      cableVolume.setVisAttributes(oddd, x_cable.visStr());

      for (int side = -1; side < 2; side += 2) {
        // Place the pipe in the stave
        staveAssembly.placeVolume(
            cableVolume,
            Transform3D(
                RotationX(0.5 * M_PI),
                Position(
                    x_cable.x_offset() + 2.05 * modCable * x_cable.rmax(),
                    side * (staveHlength - 0.5 * cableLength + x_cable.dz()),
                    x_cable.z_offset())));
      }
    }
  }
}

/// Standard create_element(...) DD4hep method for Pixel barrel
///
/// @param oddd the detector to which this is addedded
/// @param xml the input xml element
/// @param sens the sensitive detector descrition
///
/// @return a reference counted DetElement
static Ref_t create_element(Detector &oddd, xml_h xml, SensitiveDetector sens) {
  xml_det_t x_det = xml;
  string detName = x_det.nameStr();

  // Make DetElement
  DetElement barrelDetector(detName, x_det.id());
  dd4hep::xml::setDetectorTypeFlag(xml, barrelDetector);

  auto &params = ODDHelper::ensureExtension<dd4hep::rec::VariantParameters>(
      barrelDetector);

  // Add Extension to DetElement for the RecoGeometry
  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    ODDHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params,
                                         "boundary_material");
  }

  // Make Volume
  dd4hep::xml::Dimension x_det_dim(x_det.dimensions());
  string barrelShapeName = x_det_dim.nameStr();

  // The Shape and Volume
  Tube barrelShape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.dz());
  Volume barrelVolume(detName, barrelShape, oddd.air());
  barrelVolume.setVisAttributes(oddd, x_det.visStr());

  // Create the stave volume and DetElement tree
  xml_comp_t x_stave = x_det.child(_U(stave));
  Assembly staveAssembly("stave");
  // Visualization
  staveAssembly.setVisAttributes(oddd, x_stave.visStr());
  // DetElement tree
  DetElement staveElementTemplate("StaveElementTemplate", 0);

  // Create the module components
  xml_comp_t x_module = x_det.child(_U(module));
  double ylength = 0.;
  auto module =
      ODDModuleHelper::assembleRectangularModule(oddd, sens, x_module, ylength);

  // Place the modules into the stave
  double gap = x_stave.gap();
  unsigned int nModules = x_stave.nmodules();
  double ystep = ylength + gap;
  double ymin = (nModules * 0.5 - 0.5) * ylength;
  double staveHlength = ymin + 0.5 * ylength;

  // Loop over the modules and place them in the stave
  for (unsigned int moduleNum = 0; moduleNum < nModules; ++moduleNum) {
    // Place them along local y
    PlacedVolume placedModule = staveAssembly.placeVolume(
        module.first, Position(0., -ymin + moduleNum * ystep, 0.));
    placedModule.addPhysVolID("module", moduleNum);

    string moduleName = _toString((int)moduleNum, "module%d");
    // Clone the detector element
    auto moduleElement = module.second.clone(moduleName, moduleNum);
    moduleElement.setPlacement(placedModule);
    // Assign it as child to the stave template
    staveElementTemplate.add(moduleElement);
  }

  // Complete the staveStructure
  completeStaveStructure(oddd, x_stave, staveAssembly, staveHlength, ylength);

  std::vector<double> layerR;

  // Loop over the layers to build staves
  size_t layerNum = 0;
  for (xml_coll_t lay(xml, _U(layer)); lay; ++lay, ++layerNum) {
    xml_comp_t x_layer = lay;

    // The Layer envelope volume
    string layerName = detName + std::to_string(layerNum);
    Volume layerVolume(
        layerName,
        Tube(x_layer.rmin(), x_layer.rmax(), staveHlength + x_layer.outer_z()),
        oddd.air());
    // Visualization
    layerVolume.setVisAttributes(oddd, x_layer.visStr());

    // The DetElement tree, keep it flat
    DetElement layerElement(barrelDetector, x_layer.nameStr(), layerNum);
    auto &layerParams =
        ODDHelper::ensureExtension<dd4hep::rec::VariantParameters>(
            layerElement);

    // Place the staves in the layer
    unsigned int nStaves = x_layer.nphi();
    double phiStep = 2. * M_PI / nStaves;
    double phiTilt = x_layer.phi_tilt();
    double phi0 = x_layer.phi0();
    double r = x_layer.r();
    layerR.push_back(r);

    // Loop over the staves and place them
    for (unsigned int staveNum = 0; staveNum < nStaves; ++staveNum) {
      string staveName = _toString((int)staveNum, "stave%d");
      // position of the stave
      double phi = phi0 + staveNum * phiStep;
      double x = r * cos(phi);
      double y = r * sin(phi);
      // Now place the stave
      PlacedVolume placedStave = layerVolume.placeVolume(
          staveAssembly,
          Transform3D(RotationY(0.5 * M_PI) * RotationZ(0.5 * M_PI) *
                          RotationY(phi + phiTilt),
                      Position(x, y, 0.)));
      placedStave.addPhysVolID("stave", staveNum);

      // Clone the stave element from the template
      DetElement staveElement = staveElementTemplate.clone(staveName, staveNum);
      staveElement.setPlacement(placedStave);
      // Add to the layer element
      layerElement.add(staveElement);
    }

    // Place the support cylinder per layer
    buildSupportCylinder(oddd, layerVolume, x_layer, layerR);

    // Cleanup the templates
    /// @todo check if this is needed
    // dd4hep::detail::destroyHandle(staveElementTemplate);
    // dd4hep::detail::destroyHandle(module.second);

    layerParams.set<double>("envelope_r_min", 1.);
    layerParams.set<double>("envelope_r_max", 5.);
    layerParams.set<double>("envelope_z_min", 5.);
    layerParams.set<double>("envelope_z_max", 5.);
    // Add the proto layer material
    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      ODDHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams,
                                           "layer_material");
    }

    PlacedVolume placedLayer = barrelVolume.placeVolume(layerVolume);
    placedLayer.addPhysVolID("layer", layerNum);

    // Assign layer DetElement to layer volume
    layerElement.setPlacement(placedLayer);
  }

  // Place the additional support cylinders per detector
  buildSupportCylinder(oddd, barrelVolume, x_det, layerR);

  if (x_det.hasChild(_Unicode(services))) {
    // Grab the services
    xml_comp_t x_services = x_det.child(_Unicode(services));
    if (x_services.hasChild(_Unicode(cable_routing))) {
      xml_comp_t x_cable_routing = x_services.child(_Unicode(cable_routing));
      buildBarrelRouting(oddd, barrelVolume, x_cable_routing, layerR);
    }
    if (x_services.hasChild(_Unicode(cooling_routing))) {
      xml_comp_t x_cooling_routing =
          x_services.child(_Unicode(cooling_routing));
      buildBarrelRouting(oddd, barrelVolume, x_cooling_routing, layerR);
    }
  }

  // Place Volume
  Volume motherVolume = oddd.pickMotherVolume(barrelDetector);
  Position translation(0., 0., x_det_dim.z());
  PlacedVolume placedBarrel =
      motherVolume.placeVolume(barrelVolume, translation);
  // "system" is hard coded in the DD4Hep::VolumeManager
  placedBarrel.addPhysVolID("system", barrelDetector.id());
  barrelDetector.setPlacement(placedBarrel);

  // And return the detector element
  return barrelDetector;
}

DECLARE_DETELEMENT(ODDPixelBarrel, create_element)
