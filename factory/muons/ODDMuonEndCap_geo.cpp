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
		xml_comp_t x_gas = x_l.child(_Unicode(tubs));
		xml_comp_t x_shell = x_l.child(_Unicode(shell));

		if (x_gas.isSensitive()){
			sens.setType("tracker");
		}

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

			// the radial position of the chamber
			r += x_ch.dy();

			// shape of the chamber and the tubes of each layer
			Tube driftTubeShape(x_gas.rmin(), x_gas.rmax(), x_ch.dy());
			Tube shellTubeShape(x_shell.rmin(), x_shell.rmax(), x_ch.dy());
			Box chBox(x_ch.dx(), x_ch.dy(), x_ch.dz());

			for (int i = 0; i < x_l.nphi(); i++){

				double phi = phi0 + i * phistep;

				// position of chamber
				double z = (i % 2) ? x_ch.zmax() : x_ch.zmin();
				double x = r * cos(phi);
				double y = r * sin(phi);

				// create the volumes for the chambers and place them
				chamberName = _toString((int)chamberNum, "chamber%d");
				Volume chVolume(chamberName, chBox, oddd.air());
				chVolume.setVisAttributes(oddd, x_ch.visStr());
				PlacedVolume pvchamber = layerMuonVolume.placeVolume(chVolume, Transform3D(RotationZ(i * phistep), Position(x, y, z)));
				pvchamber.addPhysVolID("chamber", chamberNum++);
				DetElement chamberElement(layerElement, chamberName, chamberNum);
				chamberElement.setPlacement(pvchamber);

				// loop over the tube layers along z
				double zt = -x_ch.dz() + x_shell.rmax();
				for (int iz = 0; iz < x_gas.nz(); iz++){
					int ntubesx = x_ch.dx() / (x_shell.rmax());
					double xt = -x_ch.dx() + x_gas.x_offset();

					if (iz % 2 != 0){
						xt += x_gas.x_offset();
						ntubesx -= 1;
					}
					// loop over the tubes along x
					for (int ix = 0; ix < ntubesx; ix++){

						// create and place the tubes inside the chambers (gas+shell)
						tubeName = _toString((int)tubeNum, "tube%d");
						Volume driftTubeVolume(tubeName, driftTubeShape, oddd.air());
						driftTubeVolume.setVisAttributes(oddd, x_gas.visStr());
						driftTubeVolume.setSensitiveDetector(sens);
						PlacedVolume pvtube = chVolume.placeVolume(driftTubeVolume, Transform3D(RotationX(0.5 * M_PI), Position(xt, 0, zt)));
						pvtube.addPhysVolID("tube", tubeNum++);
						DetElement tubeElement(chamberElement, tubeName, tubeNum);
						tubeElement.setPlacement(pvtube);

						Volume shellTubeVolume(x_shell.nameStr(), shellTubeShape, oddd.air());
						shellTubeVolume.setVisAttributes(oddd, x_shell.visStr());
						chVolume.placeVolume(shellTubeVolume, Transform3D(RotationX(0.5 * M_PI), Position(xt, 0, zt)));

						// shift for the next tube -along x
						xt += 2 * x_gas.x_offset();
					}
					// shift for the next tubes layer -along z
					zt += x_gas.z_offset();
				}
			}

			// reach the top plane of the chamber -to move to the next chamber
			r += x_ch.dy();
		}
	}

	// visualize the barrel cylinder
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
