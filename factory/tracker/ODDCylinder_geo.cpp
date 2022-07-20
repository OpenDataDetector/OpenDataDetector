// Open Data Dector project
//
// (c) 2021 CERN for the benefit of the ODD project
//
// Mozilla Public License Version 2.0

#include "DD4hep/DetFactoryHelper.h"
#include "ODDHelper.hpp"
#include "XML/Utilities.h"

using namespace std;
using namespace dd4hep;

/// Standard create_element(...) create a simple cylinder
///
/// @param oddd the detector to which this is addedded
/// @param xml the input xml element
/// @param sens is ignored
///
/// @return a reference counted DetElement
static Ref_t create_element(Detector &oddd, xml_h xml,
                            SensitiveDetector /*sens*/) {
  xml_det_t x_det = xml;
  string detName = x_det.nameStr();

  // Make Volume
  xml_comp_t x_det_tubs = x_det.child(_U(tubs));

  // Make DetElement
  DetElement cylinderElement(detName, x_det.id());
  dd4hep::xml::setDetectorTypeFlag(xml, cylinderElement);

  auto &params = ODDHelper::ensureExtension<dd4hep::rec::VariantParameters>(
      cylinderElement);

  // Add the proto boundary material
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    ODDHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params,
                                         "boundary_material");
  }

  // Add the proto layer material
  for (xml_coll_t lmat(x_det_tubs, _Unicode(layer_material)); lmat; ++lmat) {
    xml_comp_t x_layer_material = lmat;
    ODDHelper::xmlToProtoSurfaceMaterial(x_layer_material, params,
                                         "layer_material");
  }

  string shapeName = x_det_tubs.nameStr();
  Tube tubeShape(shapeName, x_det_tubs.rmin(), x_det_tubs.rmax(),
                 x_det_tubs.dz());
  Volume tubeVolume(detName, tubeShape,
                    oddd.material(x_det_tubs.materialStr()));
  tubeVolume.setVisAttributes(oddd, x_det.visStr());

  // Place it in the mother
  Volume motherVolume = oddd.pickMotherVolume(cylinderElement);
  PlacedVolume placedTube = motherVolume.placeVolume(tubeVolume);
  placedTube.addPhysVolID(detName, cylinderElement.id());
  cylinderElement.setPlacement(placedTube);

  // And return the element for further parsing
  return cylinderElement;
}

DECLARE_DETELEMENT(ODDCylinder, create_element)
