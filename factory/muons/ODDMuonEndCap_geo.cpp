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
	
	xml_comp_t x_tb = x_det.child(_U(tubs));
		
	double phi0 = 0.5*M_PI;
	double phistep = 2*M_PI/x_det_dim.nphi();
	double rmin = x_det_dim.rmin();
	double r = rmin;
	double rtubemin = x_tb.rmin();
	double rtubemax = x_tb.rmax();

	//The shape and volume 
	Tube endcapMuonShape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.dz());
	Volume endcapMuonVolume(detName, endcapMuonShape,oddd.air());

	size_t chamberNum = 0;

	for(xml_coll_t chamber(x_det, _Unicode(chamber)); chamber; chamber++, chamberNum++){

		xml_comp_t x_ch = chamber;

		//the radial position of the chamber
		r += x_ch.dy();

		for(int i=0; i<x_det_dim.nphi(); i++){

			double phi = phi0 + i*phistep;

			//position of chamber on z
			double z = x_ch.zmax();

			if(i%2==0){

				 z = x_ch.zmin();				 
			}

			Box chBox(x_ch.dx(), x_ch.dy(), x_ch.dz());
			Volume chVolume(x_ch.nameStr(), chBox, oddd.air());
			chVolume.setVisAttributes(oddd,x_ch.visStr());	

			//position of the chamber volumes
			double x = r*cos(phi);
			double y = r*sin(phi);

			//place the tubes inside the chambers
			Tube driftTubeShape(rtubemin, rtubemax, x_ch.dy());
			Volume driftTubeVolume(x_tb.nameStr(), driftTubeShape, oddd.air());
			driftTubeVolume.setVisAttributes(oddd, x_tb.visStr());

			//the number of tubes along -x included the offset of one radius
			int ntubesx = (x_ch.dx()-rtubemax)/(rtubemax);
			//loop over the layers along z
			double zt = -x_ch.dz() + rtubemax;
			for(int iz=0; iz<x_tb.nz(); iz++){

				double xt = -x_ch.dx() + x_tb.x_offset();
				if(iz%2!=0){
					xt = -x_ch.dx() + 2*x_tb.x_offset();
				}
				//loop over the tubes along x
				for(int ix=0; ix<ntubesx; ix++){
				
					chVolume.placeVolume(driftTubeVolume, Transform3D(RotationX(0.5*M_PI),Position(xt,0,zt)));
					//shift for the next tube -along x
					xt+=2*x_tb.x_offset();

				}
				//shift for the next tubes layer -along z
				zt+=x_tb.z_offset();

			}

			endcapMuonVolume.placeVolume(chVolume, Transform3D(RotationZ(i*phistep),Position(x,y,z)));


		}

		//reach the top plane of the chamber -to move to the next chamber
		r+=x_ch.dy();


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
