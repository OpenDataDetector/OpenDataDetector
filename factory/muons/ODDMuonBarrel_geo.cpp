
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

	//The shape and volume 
	Tube barrelMuonShape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.dz());
	Volume barrelMuonVolume(detName, barrelMuonShape,oddd.air());
	
	xml_comp_t x_layer = x_det.child(_Unicode(layer));
    xml_comp_t x_ch = x_layer.child(_Unicode(chamber));
    xml_comp_t x_ml = x_ch.child(_Unicode(multilayer));
    xml_comp_t x_tb = x_ml.child(_U(tubs));
    xml_comp_t x_tlayer = x_ml.child(_U(layer));
    
    //Parameters for the volumes
    double rmin = x_layer.rmin();
    double rmax = sqrt(pow(0.5*x_ch.dx()+0.4*x_ch.dx(),2)+pow(rmin+x_ch.dy()+0.5*x_ch.dy(),2));
    
	
	for(int l=0; l< 3; l++){        
        
		Tube layerMuonShape(rmin, rmax, x_layer.dz());
		Volume layerVolume(x_layer.nameStr() + std::to_string(l), layerMuonShape, oddd.air());
		layerVolume.setVisAttributes(oddd, x_layer.visStr());  
        
        //The radial distance on xy for the big and small barrel chambers
		double rb = sqrt(pow(x_ch.x(),2) + pow(rmin+0.5*x_ch.dy(),2));
		double rsm = sqrt(pow(x_ch.x(),2) + pow(rmin+0.5*x_ch.dy() + x_ch.dy(),2));

        double phistep = 2*M_PI/x_layer.nphi();
        double phi0 = 0.5*M_PI;
    
        
        //place the chambers for each layer
        for(int i=0; i<x_layer.nphi(); i++){
            double phi = phi0 + i*phistep;
			
			//z position for side A and side B
			double za = -0.25*x_ch.dx();
			double zb = -za;
			double x,y;

			double zstep = x_ch.dx() + 30.;
			

			for(int j=0; j<0.5*x_layer.nz(); j++){

			za = za + zstep;
			zb = zb - zstep;

			//Big barrel chambers
			if(i%2==0){

            //chamber
			Box chBox(0.5*(x_ch.dx()+l*40),0.5*x_ch.dy(),0.5*(x_ch.dz()+l*30));
			Volume chVolume("MDT_Chamber_Big", chBox, oddd.air());
			chVolume.setVisAttributes(oddd,x_ch.visStr());	

			x = rb*cos(phi);
			y = rb*sin(phi);
            
            //place the multilayers with the tubes 
            //multilayer
            Box mlBox(chBox.x(), 0.4*chBox.y(), chBox.z());
            Volume mlVolume("MDT_MultiLayer", mlBox, oddd.air());
            mlVolume.setVisAttributes(oddd,x_ml.visStr());           
            int ntubes = mlBox.x()/x_tb.rmax();           
            Tube driftTubeShape(x_tb.rmin(), x_tb.rmax(), mlBox.z());
            Volume driftTubeVolume(x_tb.nameStr(), driftTubeShape, oddd.air());
            driftTubeVolume.setVisAttributes(oddd,x_tb.visStr());
            
             for(int nl=0; nl<x_tlayer.repeat(); nl++){                
                
                for(int t=0; t<ntubes; t++){
                    
                    mlVolume.placeVolume(driftTubeVolume,Transform3D(Position(-mlBox.x()+(2*t+1)*x_tb.rmax(), -mlBox.y()+(2*nl+1)*x_tb.rmax(),0 )));   
                  
                }
            }
            
            chVolume.placeVolume(mlVolume, Transform3D(Position(0, - chBox.y() + mlBox.y(),0))); 
            chVolume.placeVolume(mlVolume, Transform3D(Position(0, + chBox.y() - mlBox.y(),0)));

			layerVolume.placeVolume(chVolume, Transform3D(RotationZ(i*phistep),Position(x,y,za)));
			layerVolume.placeVolume(chVolume, Transform3D(RotationZ(i*phistep),Position(x,y,zb)));
		

			}else{

			Box chBox(0.4*(x_ch.dx()+l*40),0.25*x_ch.dy(),0.4*(x_ch.dz()+l*30));
			Volume chVolume("MDT_Chamber_Small", chBox, oddd.air());
			chVolume.setVisAttributes(oddd,x_ch.visStr());

			 x = rsm*cos(phi);
			 y = rsm*sin(phi);
             
             //place the multilayers with the tubes 
            //multilayer
            Box mlBox(chBox.x(),0.4*chBox.y(), chBox.z());
            Volume mlVolume("MDT_MultiLayer", mlBox, oddd.air());
            mlVolume.setVisAttributes(oddd,x_ml.visStr());
            Tube driftTubeShape(0.5*x_tb.rmin(), 0.5*x_tb.rmax(), mlBox.z());
            Volume driftTubeVolume(x_tb.nameStr(), driftTubeShape, oddd.air());
            driftTubeVolume.setVisAttributes(oddd,x_tb.visStr());
            int ntubes = 2*mlBox.x()/x_tb.rmax();
             for(int nl=0; nl<x_tlayer.repeat(); nl++){           
                
                for(int t=0; t<ntubes; t++){
                    
                    mlVolume.placeVolume(driftTubeVolume,Transform3D(Position(-mlBox.x()+(2*t+1)*0.5*x_tb.rmax(), -mlBox.y()+(2*nl+1)*0.5*x_tb.rmax(),0)));   
                  
                }
            }
             
            chVolume.placeVolume(mlVolume, Transform3D(Position(0, - chBox.y() + mlBox.y(),0))); 
            chVolume.placeVolume(mlVolume, Transform3D(Position(0, + chBox.y() - mlBox.y(),0)));


			layerVolume.placeVolume(chVolume, Transform3D(RotationZ(i*phistep),Position(x,y,za)));
			layerVolume.placeVolume(chVolume, Transform3D(RotationZ(i*phistep),Position(x,y,zb)));

			}

            
            
        }
    }       
        
        barrelMuonVolume.placeVolume(layerVolume, Position(0.,0.,0.));		
        
        rmin = rmax + 50.;
        rmax = sqrt(pow(0.5*x_ch.dx()+0.4*x_ch.dx(),2)+pow(rmin+x_ch.dy()+0.5*x_ch.dy(),2)); 
        std::cout<<rmax<<std::endl;
             
		
		
		
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


