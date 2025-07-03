#include "DD4hep/DetFactoryHelper.h"
#include "DDRec/DetectorData.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

using dd4hep::rec::LayeredCalorimeterData;

static void placeStaves(DetElement& parent, DetElement& stave, double rmin, int numsides, double total_thickness,
                        Volume envelopeVolume, double innerAngle, Volume sectVolume) {
  double innerRotation = innerAngle;
  double offsetRotation = -innerRotation / 2;
  double sectCenterRadius = rmin + total_thickness / 2;
  double rotX = M_PI / 2;
  double rotY = -offsetRotation;
  double posX = -sectCenterRadius * std::sin(rotY);
  double posY = sectCenterRadius * std::cos(rotY);

  for (int module = 1; module <= numsides; ++module) {
    DetElement det = module > 1 ? stave.clone(_toString(module,"stave%d")) : stave;
    Transform3D trafo(RotationZYX(0, rotY, rotX), Translation3D(-posX, -posY, 0));
    PlacedVolume pv = envelopeVolume.placeVolume(sectVolume,trafo);
    pv.addPhysVolID("module", module);
    det.setPlacement(pv);
    parent.add(det);
    rotY -= innerRotation;
    posX = -sectCenterRadius * std::sin(rotY);
    posY = sectCenterRadius * std::cos(rotY);
  }
}

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens) {
  xml_det_t x_det = e;
  Layering layering(x_det);
  xml_comp_t staves = x_det.staves();
  xml_dim_t dim = x_det.dimensions();
  string det_name = x_det.nameStr();
  Material air = description.air();
  double totalThickness = layering.totalThickness();
  double gap = xml_dim_t(x_det).gap();
  int numSides = dim.numsides();
  double detZ = dim.z();
  double rmin = dim.rmin();
  DetElement sdet(det_name, x_det.id());
  DetElement stave("stave1", x_det.id());
  Volume motherVol = description.pickMotherVolume(sdet);

  PolyhedraRegular polyhedra(numSides, rmin, rmin + totalThickness, detZ);
  Volume envelopeVol(det_name, polyhedra, air);

  // Add the subdetector envelope to the structure.
  double innerAngle = 2 * M_PI / numSides;
  double halfInnerAngle = innerAngle / 2;
  double tan_inner = std::tan(halfInnerAngle) * 2;
  double innerFaceLen = rmin * tan_inner;
  double outerFaceLen = (rmin + totalThickness) * tan_inner;
  double staveThickness = totalThickness;

  Trapezoid staveTrdOuter(innerFaceLen / 2, outerFaceLen / 2, detZ / 2, detZ / 2, staveThickness / 2);
  Volume staveOuterVol("stave_outer", staveTrdOuter, air);

  Trapezoid staveTrdInner(innerFaceLen / 2 - gap, outerFaceLen / 2 - gap, detZ / 2, detZ / 2, staveThickness / 2);
  Volume staveInnerVol("stave_inner", staveTrdInner, air);

  double layerOuterAngle = (M_PI - innerAngle) / 2;
  double layerInnerAngle = (M_PI / 2 - layerOuterAngle);
  double layer_pos_z = -(staveThickness / 2);
  double layer_dim_x = innerFaceLen / 2 - gap * 2;
  int layer_num = 1;

  // Copy of the LayeredCaloData from CLD factory
  // Create caloData object to extend driver with data required for reconstruction
  LayeredCalorimeterData* caloData = new LayeredCalorimeterData;
  caloData->layoutType = LayeredCalorimeterData::BarrelLayout;
  Segmentation seg = sens.readout().segmentation();
  std::vector<double> cellSizeVector =
      seg.segmentation()->cellDimensions(0); // Assume uniform cell sizes, provide dummy cellID
  double cell_sizeX = cellSizeVector[0];
  double cell_sizeY = cellSizeVector[1];
  caloData->inner_symmetry = numSides;
  caloData->outer_symmetry = numSides;
  /** NOTE: phi0=0 means lower face flat parallel to experimental floor
   *  This is achieved by rotating the modules with respect to the envelope
   *  which is assumed to be a Polyhedron and has its axes rotated with respect
   *  to the world by 180/nsides. In any other case (e.g. if you want to have
   *  a tip of the calorimeter touching the ground) this value needs to be computed
   */
  caloData->inner_phi0 = 0.;
  caloData->outer_phi0 = 0.;
  caloData->gap0 = 0.;
  caloData->gap1 = 0.;
  caloData->gap2 = 0.;
  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = rmin;
  caloData->extent[1] = rmin + totalThickness;
  caloData->extent[2] = 0.;
  caloData->extent[3] = detZ / 2.0;

  for (xml_coll_t xc(x_det, _U(layer)); xc; ++xc) {
    xml_comp_t x_layer = xc;
    int repeat = x_layer.repeat();            // Get number of times to repeat this layer.
    const Layer* lay = layering.layer(layer_num - 1); // Get the layer from the layering engine.
    LayeredCalorimeterData::Layer caloLayer;
    caloLayer.cellSize0 = cell_sizeX;
    caloLayer.cellSize1 = cell_sizeY;

    // Loop over repeats for this layer.
    for (int j = 0; j < repeat; j++) {
      string layer_name = _toString(layer_num, "layer%d");
      double layer_thickness = lay->thickness();
      DetElement layer(stave, layer_name, layer_num);

      // Layer position in Z within the stave.
      layer_pos_z += layer_thickness / 2;
      // Layer box & volume
      Volume layer_vol(layer_name, Box(layer_dim_x, detZ / 2, layer_thickness / 2), air);

      // Create the slices (sublayers) within the layer.
      double slice_pos_z = -(layer_thickness / 2);
      int slice_number = 1;
      double nRadiationLengths = 0.;
      double nInteractionLengths = 0.;
      double thickness_sum = 0;

      for (xml_coll_t xk(x_layer, _U(slice)); xk; ++xk) {
        xml_comp_t x_slice = xk;
        string slice_name = _toString(slice_number, "slice%d");
        double slice_thickness = x_slice.thickness();
        Material slice_material = description.material(x_slice.materialStr());
        DetElement slice(layer, slice_name, slice_number);

        slice_pos_z += slice_thickness / 2;
        // Slice volume & box
        Volume slice_vol(slice_name, Box(layer_dim_x, detZ / 2, slice_thickness / 2), slice_material);

        nRadiationLengths += slice_thickness / (2. * slice_material.radLength());
        nInteractionLengths += slice_thickness / (2. * slice_material.intLength());
        thickness_sum += slice_thickness / 2;
        // Store "inner" quantities
        caloLayer.inner_nRadiationLengths = nRadiationLengths;
        caloLayer.inner_nInteractionLengths = nInteractionLengths;
        caloLayer.inner_thickness = thickness_sum;
        caloLayer.sensitive_thickness = slice_thickness;
        nRadiationLengths += slice_thickness / (2. * slice_material.radLength());
        nInteractionLengths += slice_thickness / (2. * slice_material.intLength());
        thickness_sum += slice_thickness / 2;

        if (x_slice.isSensitive()) {
          sens.setType("calorimeter");
          slice_vol.setSensitiveDetector(sens);
        }
        // Set region, limitset, and vis.
        slice_vol.setAttributes(description, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
        // slice PlacedVolume
        PlacedVolume slice_phv = layer_vol.placeVolume(slice_vol, Position(0, 0, slice_pos_z));
        slice_phv.addPhysVolID("slice", slice_number);

        slice.setPlacement(slice_phv);
        // Increment Z position for next slice.
        slice_pos_z += slice_thickness / 2;
        // Increment slice number.
        ++slice_number;
      }
      // Store "outer" quantities
      caloLayer.outer_nRadiationLengths = nRadiationLengths;
      caloLayer.outer_nInteractionLengths = nInteractionLengths;
      caloLayer.outer_thickness = thickness_sum;

      // Set region, limitset, and vis.
      layer_vol.setAttributes(description, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());

      // Layer physical volume.
      PlacedVolume layer_phv = staveInnerVol.placeVolume(layer_vol, Position(0, 0, layer_pos_z));
      layer_phv.addPhysVolID("layer", layer_num);
      layer.setPlacement(layer_phv);

      // The rest of the data is constant; only the distance needs to be updated
      // Store the position up to the inner face of the layer
      caloLayer.distance = rmin + layer_pos_z + staveThickness / 2 - layer_thickness / 2;
      std::cout << "Layer: " << layer_num << " Rmin: " << rmin << " layer_pos_z: " << layer_pos_z
                << " Dist: " << caloLayer.distance
                << " inner_thickness: " << caloLayer.inner_thickness
                << " outer_thickness: " << caloLayer.outer_thickness
                << std::endl;
      // Push back a copy to the caloData structure
      caloData->layers.push_back(caloLayer);

      // Increment the layer X dimension.
      layer_dim_x += layer_thickness * std::tan(layerInnerAngle);    // * 2;
      // Increment the layer Z position.
      layer_pos_z += layer_thickness / 2;
      // Increment the layer number.
      ++layer_num;
    }
  }

  // Add stave inner physical volume to outer stave volume.
  staveOuterVol.placeVolume(staveInnerVol);
  if ( staves )  {
    // Set the vis attributes of the outer stave section.
    stave.setVisAttributes(description, staves.visStr(), staveInnerVol);
    stave.setVisAttributes(description, staves.visStr(), staveOuterVol);
  }
  // Place the staves.
  placeStaves(sdet, stave, rmin, numSides, totalThickness, envelopeVol, innerAngle, staveOuterVol);
  // Set envelope volume attributes.
  envelopeVol.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  xml::setDetectorTypeFlag(e, sdet);

  double z_offset = dim.hasAttr(_U(z_offset)) ? dim.z_offset() : 0.0;
  Transform3D transform(RotationZ(M_PI / numSides), Translation3D(0, 0, z_offset));
  PlacedVolume env_phv = motherVol.placeVolume(envelopeVol, transform);
  env_phv.addPhysVolID("system", sdet.id());
  env_phv.addPhysVolID("barrel", 0);
  sdet.setPlacement(env_phv);

  sdet.addExtension<LayeredCalorimeterData>(caloData);
  return sdet;
}

DECLARE_DETELEMENT(ODDPolyhedraBarrelCalorimeter, create_detector)
