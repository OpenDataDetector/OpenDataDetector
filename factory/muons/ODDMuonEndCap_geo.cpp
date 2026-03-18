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

static Ref_t create_element(Detector &oddd, xml_h xml, SensitiveDetector sens){

	xml_det_t x_det = xml;
	string detName = x_det.nameStr();

	DetElement endcapMuonDetector(x_det.nameStr(), x_det.id());
	dd4hep::xml::setDetectorTypeFlag(xml, endcapMuonDetector);

	// Make Volume
	dd4hep::xml::Dimension x_det_dim(x_det.dimensions());
	string endcapShapeName = x_det_dim.nameStr();

	// The shape and volume
	Tube endcapMuonShape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.dz());
	Volume endcapMuonVolume(detName, endcapMuonShape, oddd.air());

	unsigned int chamberNum = 1;
	unsigned int tubeNum = 1;
	unsigned int layerNum = 1;

	std::string tubeName, chamberName, layerName;

	for (xml_coll_t layer(x_det, _Unicode(layer)); layer; layer++){

		xml_comp_t x_l = layer;
	
		double phi0 = 0.5 * M_PI;
		double phistep = 2 * M_PI / x_l.nphi();
		double r = x_l.rmin();
		layerName = _toString((int)layerNum, "layer%d");
		Tube layerMuonShape(x_l.rmin(), x_l.rmax(), x_l.dz());
		Volume layerMuonVolume(layerName, layerMuonShape, oddd.air());
		PlacedVolume pvlayer = endcapMuonVolume.placeVolume(layerMuonVolume, Position(0, 0, x_l.z()));
		pvlayer.addPhysVolID("layer", layerNum++);
		DetElement layerElement(endcapMuonDetector, layerName, layerNum);
		layerElement.setPlacement(pvlayer);

		for (xml_coll_t chamber(x_l, _Unicode(chamber)); chamber; chamber++){

			xml_comp_t x_ch = chamber;
			xml_comp_t x_gas = x_ch.child(_Unicode(tubs));
			xml_comp_t x_shell = x_ch.child(_Unicode(shell));

			if (x_gas.isSensitive()){
				sens.setType("tracker");
			}


			// the radial position of the chamber
			r += x_ch.dy();
			double dx1 = 0.4* x_ch.dx();          // bottom half-width
			double dx2 = x_ch.dx();    // top half-width (narrower)
			double dy = x_ch.dy();
			double dz  = x_ch.dz();
			Trapezoid chTrap(dx1, dx2, dz, dz, dy);	

			for (int i = 0; i < x_l.nphi(); i++){

				double phi = phi0 + i * phistep;

				// position of chamber
				double z = (i % 2) ? x_ch.zmax() : x_ch.zmin();
				double x = r * cos(phi);
				double y = r * sin(phi);

				// create the volumes for the chambers and place them
				chamberName = _toString((int)chamberNum, "chamber%d");
				Volume chVolume(chamberName, chTrap, oddd.air());
				chVolume.setVisAttributes(oddd, x_ch.visStr());
	
				PlacedVolume pvchamber = layerMuonVolume.placeVolume(chVolume, Transform3D(RotationZ(i * phistep), Position(x, y, z))*Transform3D(RotationX(-0.5*M_PI)));
				pvchamber.addPhysVolID("chamber", chamberNum++);
				DetElement chamberElement(layerElement, chamberName, chamberNum);
				chamberElement.setPlacement(pvchamber);

				// loop over the tube layers along z
				double zt = -x_ch.dy() + x_shell.rmax();
				int ntubesy = x_ch.dz()/(x_shell.rmax()); // global y direction of the tube layers
				int ntubesz = x_ch.dy() / (x_shell.rmax()); // global z direction of the tube layers
				double lowerTubeLength = dx1;
				double upperTubeLength = dx2;
				double h{0.};
			
				for (int iz = 0; iz < ntubesz; iz++){
					double yt = -x_ch.dz() + x_shell.rmax();
					
					for (int iy = 0; iy < ntubesy; iy++){			
																
						double tubeLength = lowerTubeLength + (upperTubeLength-lowerTubeLength)*(h/(2*dy));
						
						Tube driftTubeShape(x_gas.rmin(), x_gas.rmax(), tubeLength);
						Tube shellTubeShape(x_shell.rmin(), x_shell.rmax(), tubeLength);
						// loop over the tubes along x
							
							// create and place the tubes inside the chambers (gas+shell)
						tubeName = _toString((int)tubeNum, "tube%d");
						Volume driftTubeVolume(tubeName, driftTubeShape, oddd.material(x_gas.materialStr()));
						driftTubeVolume.setVisAttributes(oddd, x_gas.visStr());
						driftTubeVolume.setSensitiveDetector(sens);
						PlacedVolume pvtube = chVolume.placeVolume(driftTubeVolume, Transform3D(RotationY(0.5 * M_PI), Position(0, yt, zt)));
						pvtube.addPhysVolID("tube", tubeNum++);
						DetElement tubeElement(chamberElement, tubeName, tubeNum);
						tubeElement.setPlacement(pvtube);

						Volume shellTubeVolume(x_shell.nameStr(), shellTubeShape, oddd.material(x_shell.materialStr()));
						shellTubeVolume.setVisAttributes(oddd, x_shell.visStr());
						chVolume.placeVolume(shellTubeVolume, Transform3D(RotationY(0.5 * M_PI), Position(0, yt, zt)));

							// shift for the next tube -along z
						yt += 2 *x_shell.rmax();
					}
						// shift for the next tubes layer -along y
					zt += 2*x_shell.rmax();
					h += 2* x_shell.rmax();
				}
			}

			// reach the top plane of the chamber -to move to the next chamber
			r += x_ch.dy();
		}
	}

	// visualize the endcap cylinder
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
