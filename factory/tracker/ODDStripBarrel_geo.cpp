// Open Data Dector project
//
// (c) 2021 CERN for the benefit of the ODD project
//
// Mozilla Public License Version 2.0

#include <vector>

#include "DD4hep/DetFactoryHelper.h"
#include "ODDHelper.hpp"
#include "ODDModuleHelper.hpp"
#include "ODDServiceHelper.hpp"
#include "XML/Utilities.h"

using namespace std;
using namespace dd4hep;

/// Standard create_element(...) DD4hep method for Strip barrel
///
/// This can build long/short strip barrels for the Open Data Detector,
/// allowing a variable number of cylindrical layers.
///
/// It supports double (with stereo angle) and single sided modules,
/// the latter are assumed to be rectangular.
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

  // Add Extension to DetElement for the RecoGeometry
  auto &params = ODDHelper::ensureExtension<dd4hep::rec::VariantParameters>(
      barrelDetector);
  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    ODDHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params,
                                         "boundary_material");
  }

  // Make Volume
  dd4hep::xml::Dimension x_det_dim(x_det.dimensions());
  string barrelShapeName = x_det_dim.nameStr();

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

  // Build a template module for the Barrel
  xml_comp_t x_module = x_det.child(_U(module));
  double length = 0.;
  auto module =
      ODDModuleHelper::assembleRectangularModule(oddd, sens, x_module, length);

  // Place the modules into the stave
  double gap = x_stave.gap();
  unsigned int nModules = x_stave.nmodules();
  double ystep = length + gap;
  double ymin = (nModules * 0.5 - 0.5) * length;
  double staveHlength = ymin + 0.5 * length;

  // Loop over the modules and place them in the stave
  for (unsigned int moduleNum = 0; moduleNum < nModules; ++moduleNum) {
    double positionY = -ymin + moduleNum * ystep;
    // Place the cable bundle, one per stave
    if (x_stave.hasChild(_U(eltube))) {
      // Retrieve cable parameters
      xml_comp_t x_cable = x_stave.child(_U(eltube));

      double rMin = x_cable.rmin();
      double rMax = x_cable.rmax();

      // For an odd number of modules this will create an asymmetric powering
      // (as it should)
      double rStep = (rMax - rMin) / (0.5 * nModules);
      double rCable = rMin + abs(moduleNum - 0.5 * nModules) * rStep;

      Tube cable(0., rCable, 0.495 * ystep);
      // Create the scable volume
      Volume cableVolume("Cable", cable, oddd.material(x_cable.materialStr()));
      cableVolume.setVisAttributes(oddd, x_cable.visStr());

      // Place the pipe in the stave
      staveAssembly.placeVolume(
          cableVolume, Transform3D(RotationX(0.5 * M_PI),
                                   Position(x_cable.x_offset(), positionY,
                                            x_cable.z_offset())));
    }

    // Place them along local y
    PlacedVolume placedModule =
        staveAssembly.placeVolume(module.first, Position(0., positionY, 0.));
    placedModule.addPhysVolID("module", moduleNum);

    string moduleName = _toString((int)moduleNum, "module%d");
    // Clone the detector element
    auto moduleElement = module.second.clone(moduleName, moduleNum);
    moduleElement.setPlacement(placedModule);
    // Assign it as child to the stave template
    staveElementTemplate.add(moduleElement);
  }

  // Remember the layer radii
  std::vector<double> layerR;

  // Loop over the layers to build staves
  size_t layerNum = 0;
  for (xml_coll_t lay(xml, _U(layer)); lay; ++lay, ++layerNum) {
    xml_comp_t x_layer = lay;

    string layerName = detName + std::to_string(layerNum);
    // The Module envelope volume
    Volume layerVolume(
        layerName,
        Tube(x_layer.rmin(), x_layer.rmax(), staveHlength + x_layer.outer_z()),
        oddd.air());
    // Visualization
    layerVolume.setVisAttributes(oddd, x_layer.visStr());

    // The DetElement tree, keep it flat
    DetElement layerElement(barrelDetector, x_layer.nameStr(), layerNum);

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

    // Place the support cylinder
    std::vector<double> dummyR;
    buildSupportCylinder(oddd, barrelVolume, x_layer, dummyR);

    auto &layerParams =
        ODDHelper::ensureExtension<dd4hep::rec::VariantParameters>(
            layerElement);

    layerParams.set<double>("envelope_r_min", 10.);
    layerParams.set<double>("envelope_r_max", 25.);
    layerParams.set<double>("envelope_z_min", 10.);
    layerParams.set<double>("envelope_z_max", 10.);

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

  }  // loop over layers

  // Place the support rails
  buildSupportCylinder(oddd, barrelVolume, x_det, layerR);

  // Route the services out on both sides
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

  // And return it
  return barrelDetector;
}

DECLARE_DETELEMENT(ODDStripBarrel, create_element)
