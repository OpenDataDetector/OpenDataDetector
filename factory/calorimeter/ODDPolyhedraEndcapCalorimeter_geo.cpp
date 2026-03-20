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

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {
  xml_det_t   x_det     = e;
  xml_comp_t  staves    = x_det.staves();
  xml_dim_t   dim       = x_det.dimensions();
  int         det_id    = x_det.id();
  bool        reflect   = x_det.reflect(true);
  string      det_name  = x_det.nameStr();
  Material    air       = description.air();
  double      gap       = xml_dim_t(x_det).gap();
  int         numsides  = dim.numsides();
  double      rmin      = dim.rmin();
  double      rmax      = dim.rmax()*std::cos(M_PI/numsides);
  double      zmin      = dim.zmin();
  Layering    layering(x_det);
  double      totalThickness = layering.totalThickness();
  double      detZ      = totalThickness;
  Volume      endcapVol(det_name+"_endcap",PolyhedraRegular(numsides,rmin,rmax,totalThickness),air);
  DetElement  endcap(det_name+"_endcap",det_id);
  DetElement  stave("stave1", x_det.id());

  // Place modules/staves inside the polyhedra, as for the barrel detector
  double innerAngle = 2 * M_PI / numsides;
  double halfInnerAngle = innerAngle / 2;
  double tan_inner = std::tan(halfInnerAngle) * 2;
  double innerFaceLen = rmin * tan_inner;
  double staveThickness = rmax - rmin;
  double outerFaceLen = (rmin + staveThickness) * tan_inner;

  Trapezoid staveTrdOuter(innerFaceLen / 2, outerFaceLen / 2, detZ / 2, detZ / 2, staveThickness / 2);
  Volume staveOuterVol("stave_outer", staveTrdOuter, air);

  Trapezoid staveTrdInner(innerFaceLen / 2 - gap, outerFaceLen / 2 - gap, detZ / 2, detZ / 2, staveThickness / 2);
  Volume staveInnerVol("stave_inner", staveTrdInner, air);

  double layer_pos_y = -(detZ / 2);
  size_t layer_num = 1;

  // Copy of the LayeredCaloData from CLD factory
  // Create caloData object to extend driver with data required for reconstruction
  LayeredCalorimeterData* caloData = new LayeredCalorimeterData;
  caloData->layoutType = LayeredCalorimeterData::BarrelLayout;
  Segmentation seg = sens.readout().segmentation();
  std::vector<double> cellSizeVector =
      seg.segmentation()->cellDimensions(0); // Assume uniform cell sizes, provide dummy cellID
  double cell_sizeX = cellSizeVector[0];
  double cell_sizeY = cellSizeVector[1];
  caloData->inner_symmetry = numsides;
  caloData->outer_symmetry = numsides;
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

  endcapVol.setAttributes(description,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());

  // difference from the barrel: layers are placed along z, not along radius
  for(xml_coll_t xc(x_det,_U(layer)); xc; ++xc)  {
    xml_comp_t       x_layer  = xc;
    int              l_repeat = x_layer.repeat();
    std::cout << "ODD Polyhedra ENDCAP.\n Number of layers: " << l_repeat << std::endl;
    if ( l_repeat <= 0 ) throw std::runtime_error(x_det.nameStr()+"> Invalid repeat value");
    const Layer* lay = layering.layer(layer_num - 1); // Get the layer from the layering engine.
    LayeredCalorimeterData::Layer caloLayer;
    caloLayer.cellSize0 = cell_sizeX;
    caloLayer.cellSize1 = cell_sizeY;

    for (int j = 0; j < l_repeat; j++) {
      string layer_name = _toString(static_cast<int>(layer_num), "layer%d");
      std::cout << "- layer named " << layer_name << std::endl;
      double layer_thickness = lay->thickness();
      DetElement layer(stave, layer_name, static_cast<int>(layer_num));
      layer_pos_y += layer_thickness / 2;
      // Layer trapezoid shape & volume
      Trapezoid layer_trd(innerFaceLen / 2 - gap, outerFaceLen / 2 - gap, layer_thickness / 2, layer_thickness / 2, staveThickness / 2);
      Volume layer_vol(layer_name, layer_trd, air);
      std::cout << "Layer of trapezoidal shape with dimensions " <<  layer_thickness / 2 << " , " << staveThickness / 2 << std::endl;

      // Create the slices (sublayers) within the layer.
      double slice_pos_y = -(layer_thickness / 2);
      int slice_number = 1;
      double nRadiationLengths = 0.;
      double nInteractionLengths = 0.;
      double thickness_sum = 0;

      for(xml_coll_t xs(x_layer,_U(slice)); xs; ++xs)  {
        xml_comp_t x_slice = xs;
        string slice_name = _toString(slice_number, "slice%d");
        double slice_thickness = x_slice.thickness();
        Material slice_material = description.material(x_slice.materialStr());
        DetElement slice(layer, slice_name, slice_number);
        slice_pos_y += slice_thickness / 2;
        // Slice volume
        Trapezoid slice_trd(innerFaceLen / 2 - gap, outerFaceLen / 2 - gap, slice_thickness/ 2, slice_thickness / 2, staveThickness / 2);
        Volume slice_vol(slice_name, slice_trd, slice_material);

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
        PlacedVolume slice_phv = layer_vol.placeVolume(slice_vol, Position(0, slice_pos_y, 0));
        slice_phv.addPhysVolID("slice", slice_number);
        std::cout << "Slice " << slice_number << " made of " << x_slice.materialStr() << " with name " << slice_name
                  << " and half sizes " << slice_thickness / 2 << " , " << staveThickness / 2 <<
                  "     placed in layer at y = " << slice_pos_y  << std::endl;

        slice.setPlacement(slice_phv);
        // Increment Z position for next slice.
        slice_pos_y += slice_thickness / 2;
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
      PlacedVolume layer_phv = staveInnerVol.placeVolume(layer_vol, Position(0, layer_pos_y, 0));
      std::cout << " Placing layer " << static_cast<int>(layer_num) << " with half dimension of " << layer_thickness / 2 
      << ", " << staveThickness / 2 << " at " << 0 << ", " << layer_pos_y <<  ", " << 0 << std::endl;
      layer_phv.addPhysVolID("layer", static_cast<int>(layer_num));
      layer.setPlacement(layer_phv);

      // The rest of the data is constant; only the distance needs to be updated
      // Store the position up to the inner face of the layer
      caloLayer.distance = rmin + layer_pos_y + staveThickness / 2 - layer_thickness / 2;
      std::cout << "Layer: " << static_cast<int>(layer_num) << " Rmin: " << rmin << " layer_pos_z: " << layer_pos_y
                << " Dist: " << caloLayer.distance
                << " inner_thickness: " << caloLayer.inner_thickness
                << " outer_thickness: " << caloLayer.outer_thickness
                << std::endl;
      // Push back a copy to the caloData structure
      caloData->layers.push_back(caloLayer);

      // Increment the layer Z position.
      layer_pos_y += layer_thickness / 2;
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
  placeStaves(endcap, stave, rmin, numsides, staveThickness, endcapVol, innerAngle, staveOuterVol);

  double z_pos = zmin+totalThickness/2;
  PlacedVolume pv;
  // Reflect it.
  if ( reflect )  {
    Assembly    assembly(det_name);
    DetElement  both_endcaps(det_name,det_id);
    Volume      motherVol = description.pickMotherVolume(both_endcaps);
    DetElement  sdetA = endcap;
    Ref_t(sdetA)->SetName((det_name+"_A").c_str());
    DetElement  sdetB = endcap.clone(det_name+"_B",x_det.id());
    xml::setDetectorTypeFlag(e, both_endcaps);
    both_endcaps.addExtension<LayeredCalorimeterData>(caloData);

    pv = assembly.placeVolume(endcapVol,Transform3D(RotationZYX(M_PI/numsides,0,0),
                                                    Position(0,0,z_pos)));
    pv.addPhysVolID("barrel", 1);
    sdetA.setPlacement(pv);

    pv = assembly.placeVolume(endcapVol,Transform3D(RotationZYX(M_PI/numsides,M_PI,0),
                                                    Position(0,0,-z_pos)));
    pv.addPhysVolID("barrel", 2);
    sdetB.setPlacement(pv);

    pv = motherVol.placeVolume(assembly);
    pv.addPhysVolID("system", det_id);
    both_endcaps.setPlacement(pv);
    both_endcaps.add(sdetA);
    both_endcaps.add(sdetB);
    return both_endcaps;
  }
  // else
  Volume motherVol = description.pickMotherVolume(endcap);
  xml::setDetectorTypeFlag(e, endcap);
  pv = motherVol.placeVolume(endcapVol,Transform3D(RotationZYX(M_PI/numsides,0,0),
                                                   Position(0,0,z_pos)));
  pv.addPhysVolID("system", det_id);
  pv.addPhysVolID("barrel", 1);
  endcap.setPlacement(pv);
  Ref_t(endcap)->SetName(det_name.c_str());
  endcap.addExtension<LayeredCalorimeterData>(caloData);
  return endcap;
}

DECLARE_DETELEMENT(ODDPolyhedraEndcapCalorimeter,create_detector)

