
// Open Data Dector project
//
// (c) 2021 CERN for the benefit of the ODD project
//
// Mozilla Public License Version 2.0



#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"

using namespace std;
using namespace dd4hep;

/// Standard create_element(...) create muon spectrometer barrel like geometry
///
/// @param oddd the detector to which this is addedded
/// @param xml the input xml element
/// @param sens is ignored
///
/// @return a reference counted DetElement

static Ref_t create_element(Detector &oddd, xml_h xml,  SensitiveDetector sens){
	xml_det_t x_det = xml;
	string detName = x_det.nameStr();

	// Construct the detector element 
	DetElement barrelMuonDetector(detName, x_det.id());
	dd4hep::xml::setDetectorTypeFlag(xml, barrelMuonDetector);

	// Make Volume
  	dd4hep::xml::Dimension x_det_dim(x_det.dimensions());
  	string barrelShapeName = x_det_dim.nameStr();	
	
    xml_comp_t x_ch = x_det.child(_Unicode(chamber));
    xml_comp_t x_ml = x_ch.child(_Unicode(multilayer));
    xml_comp_t x_tb = x_ml.child(_U(tubs));
    xml_comp_t x_tlayer = x_ml.child(_U(layer)); 

      //Parameters for the volumes
    double rmin = x_det_dim.rmin();
    //The shape and volume 
	Tube barrelMuonShape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.dz());
	Volume barrelMuonVolume(detName, barrelMuonShape,oddd.air());	   
        		
    //The radial distance on xy for the big and small barrel chambers
    //small chambers are placed above the big ones
	double rb = sqrt(pow(x_ch.x(),2) + pow(rmin+0.5*x_ch.dy(),2));
	double rsm = sqrt(pow(x_ch.x(),2) + pow(rmin+0.5*x_ch.dy() + x_ch.dy(),2));

    double phistep = 2*M_PI/x_det_dim.nphi();
    double phi0 = 0.5*M_PI; 

       
    for(int i=0; i<x_det_dim.nphi(); i++){

        double phi = phi0 + i*phistep;

        //position of the big barrel chambers

		double x=rb*cos(phi);
		double y=rb*sin(phi);

		 //z position for side A and side B
		double za = -x_ch.dz();
		double zb = -za;	

		//big  chambers dimensions
		double ch_dx = 0.5*x_ch.dx();
		double ch_dy = 0.5*x_ch.dy();
		double ch_dz = 0.5*x_ch.dz(); 
		double rtubemax = x_tb.rmax();
		double rtubemin = x_tb.rmin();
		double zstep = x_ch.dz() + x_det_dim.z_offset();
		string name = "MDT_Chamber_Big";

		//Switch to small barrel chambers when i is even
        if(i%2 == 0){

        	//position of the small barrel chambers
			x = rsm*cos(phi);
			y = rsm*sin(phi);
			//small chambers dimensions
			ch_dx*=0.8;
			ch_dy*=0.8;
			ch_dz*=0.8;
			rtubemax*=0.5;
			rtubemin*=0.5;
			name = "MDT_Chamber_Small";
		}

		for(int j=0; j<x_det_dim.nz(); j++){
			        

			//Build the chamber, the multilayer and the tube layers with the tubes
			Box chBox(ch_dx,ch_dy,ch_dz);
			Volume chVolume(name, chBox, oddd.air());
			chVolume.setVisAttributes(oddd,x_ch.visStr());	
			//Each multilayer has 0.4 of the height of the chamber
			//we need to multilayers with some space between them
			Box mlBox(ch_dx, 0.4*ch_dy, ch_dz);
            Volume mlVolume("MDT_MultiLayer", mlBox, oddd.air());
            mlVolume.setVisAttributes(oddd,x_ml.visStr());           
            int ntubes = mlBox.x()/rtubemax;
            Tube driftTubeShape(rtubemin, rtubemax, mlBox.z());
            Volume driftTubeVolume(x_tb.nameStr(), driftTubeShape, oddd.air());
            driftTubeVolume.setVisAttributes(oddd,x_tb.visStr());
            
            for(int nl=0; nl<x_tlayer.repeat(); nl++){                
                
                for(int t=0; t<ntubes; t++){

                	mlVolume.placeVolume(driftTubeVolume,Transform3D(Position(-mlBox.x()+(2*t+1)*rtubemax, -mlBox.y()+(2*nl+1)*rtubemax,0 )));   
                  
                }
            }

            chVolume.placeVolume(mlVolume, Transform3D(Position(0, - chBox.y() + mlBox.y(),0))); 
            chVolume.placeVolume(mlVolume, Transform3D(Position(0, + chBox.y() - mlBox.y(),0)));


			barrelMuonVolume.placeVolume(chVolume, Transform3D(RotationZ(i*phistep),Position(x,y,za)));
			barrelMuonVolume.placeVolume(chVolume, Transform3D(RotationZ(i*phistep),Position(x,y,zb)));

			//move along z for the next chambers
			za = za - zstep;
			zb = zb + zstep;            
            
        }

       
    }             
        

  //visualize the barrel cylinder
  barrelMuonVolume.setVisAttributes(oddd, x_det.visStr());

  // Place Volume
  Volume motherVolume = oddd.pickMotherVolume(barrelMuonDetector);
  Position translation(0., 0., x_det_dim.z());

  PlacedVolume placedMuonBarrel = motherVolume.placeVolume(barrelMuonVolume, translation);
  placedMuonBarrel.addPhysVolID("system", barrelMuonDetector.id());
  barrelMuonDetector.setPlacement(placedMuonBarrel);

  return barrelMuonDetector;
		
}

DECLARE_DETELEMENT(ODDMuonBarrel, create_element)


