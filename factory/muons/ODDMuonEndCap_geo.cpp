// Open Data Dector project
//
// (c) 2021 CERN for the benefit of the ODD project
//
// Mozilla Public License Version 2.0



#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"

using namespace std;
using namespace dd4hep;

/// Standard create_element(...) create muon spectrometer endcap like geometry
///
/// @param oddd the detector to which this is addedded
/// @param xml the input xml element
/// @param sens is ignored
///
/// @return a reference counted DetElement

static Ref_t create_element(Detector &oddd, xml_h xml,  SensitiveDetector sens){

	xml_det_t x_det = xml;
	string detName = x_det.nameStr();

	DetElement endcapMuonDetector(x_det.nameStr(), x_det.id());
	dd4hep::xml::setDetectorTypeFlag(xml, endcapMuonDetector);

	// Make Volume
  	dd4hep::xml::Dimension x_det_dim(x_det.dimensions());
  	string endcapShapeName = x_det_dim.nameStr();

	
	xml_comp_t x_layer=x_det.child(_U(layer));
	xml_comp_t x_trd = x_layer.child(_U(trd));
		
	double phi0 = 0.5*M_PI;
	double phistep = 2*M_PI/x_layer.nphi();
	double rmin = x_det_dim.rmin();

	//The shape and volume 
	Tube endcapMuonShape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.dz());
	Volume endcapMuonVolume(detName, endcapMuonShape,oddd.air());


	for(int i=0; i<x_layer.nphi(); i++){
			
			double phi = phi0 + i*phistep;

			double dz = x_trd.z1();
			double r = rmin + dz;
			double z = x_trd.zmax();
			string name="MDT_Chamber_Big";


			if(i%2==0){

				 dz = x_trd.z2();
				 r = rmin + dz;
				 z = x_trd.zmin();
				 name="MDT_Chamber_Small";
			}
			
			Trapezoid chamberShape(x_trd.x1(), x_trd.x2(), x_trd.thickness(), x_trd.thickness(), dz);
			Volume chVolume(name, chamberShape, oddd.air());
			chVolume.setVisAttributes(oddd, x_trd.visStr());
				
			double x = r*cos(phi);
			double y = r*sin(phi);

			endcapMuonVolume.placeVolume(chVolume, Transform3D(RotationZ(i*phistep)*RotationX(-0.5*M_PI),Position(x,y,z)));			
			

		}


	//visualize the barrel cylinder
  endcapMuonVolume.setVisAttributes(oddd, x_det.visStr());

  // Place Volume
  Volume motherVolume = oddd.pickMotherVolume(endcapMuonDetector);
  Position translation(0., 0., x_det_dim.z());

  PlacedVolume placedMuonEndCap = motherVolume.placeVolume(endcapMuonVolume, translation);
  placedMuonEndCap.addPhysVolID("system", endcapMuonDetector.id());
  endcapMuonDetector.setPlacement(placedMuonEndCap);

  return endcapMuonDetector;
		


}

DECLARE_DETELEMENT(ODDMuonEndCap, create_element)
