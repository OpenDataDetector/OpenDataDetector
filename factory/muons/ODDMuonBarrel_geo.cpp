
// Open Data Dector project
//
// (c) 2021 CERN for the benefit of the ODD project
//
// Mozilla Public License Version 2.0

#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"

#include <iostream>
#include <fstream>

using namespace std;
using namespace dd4hep;

/// Standard create_element(...) create muon spectrometer barrel like geometry
///
/// @param oddd the detector to which this is addedded
/// @param xml the input xml element
/// @param sens the sensitive detector
///
/// @return a reference counted DetElement

static Ref_t create_element(Detector &oddd, xml_h xml, SensitiveDetector sens)
{
	xml_det_t x_det = xml;
	string detName = x_det.nameStr();

	// Construct the detector element
	DetElement barrelMuonDetector(detName, x_det.id());
	dd4hep::xml::setDetectorTypeFlag(xml, barrelMuonDetector);

	dd4hep::xml::Dimension x_det_dim(x_det.dimensions());
	string barrelShapeName = x_det_dim.nameStr();

	// The shape and volume
	Tube barrelMuonShape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.dz());
	Volume barrelMuonVolume(detName, barrelMuonShape, oddd.air());

	// The identifiers for the chambers and the gas-tubes
	unsigned int chamberNum = 1;
	unsigned int tubeNum = 1;

	std::string tubeName, chamberName;

	for (xml_coll_t layer(x_det, _U(layer)); layer; ++layer)
	{

		xml_comp_t x_l = layer;

		xml_comp_t x_ch = x_l.child(_Unicode(chamber));
		xml_comp_t x_gas = x_ch.child(_U(tubs));
		xml_comp_t x_shell = x_ch.child(_Unicode(shell));

		if (x_gas.isSensitive())
		{
			sens.setType("tracker");
		}
		// small chambers are placed above the big ones along y
		double rb = sqrt(pow(x_ch.x(), 2) + pow(x_l.rmin() + x_ch.dy(), 2));
		double rsm = sqrt(pow(x_ch.x(), 2) + pow(x_l.rmin() + x_ch.dy() + 2 * x_ch.dy() + x_ch.y_offset(), 2));

		double phistep = 2 * M_PI / x_l.nphi();
		double phi0 = 0.5 * M_PI;

		for (int i = 0; i < x_l.nphi(); i++)
		{

			double phi = phi0 + i * phistep;

			// Define the volumes' parameters

			// position of the big barrel chambers
			double x = rb * cos(phi);
			double y = rb * sin(phi);

			// z position for side A and side B
			double za = -x_ch.dz() - x_l.z_offset();
			double zb = -za;

			// big  chambers dimensions
			double ch_dx = x_ch.dx();
			double ch_dy = x_ch.dy();
			double ch_dz = x_ch.dz();
			double rtubemax = x_gas.rmax();
			double rtubemin = x_gas.rmin();
			double rshellmin = x_shell.rmin();
			double rshellmax = x_shell.rmax();
			double zstep = 2 * x_ch.dz() + x_l.z_offset();

			// Switch to small barrel chambers when i is even
			if (i % 2 == 0)
			{

				// position of the small barrel chambers
				x = rsm * cos(phi);
				y = rsm * sin(phi);
				// small chambers dimensions
				ch_dx *= 0.9;
				ch_dy *= 0.9;
				ch_dz *= 0.9;
				rtubemax *= 0.9;
				rtubemin *= 0.9;
				rshellmin *= 0.9;
				rshellmax *= 0.9;
			}

			// Build the chamber,  and the tube layers with the tubes
			Box chBox(ch_dx, ch_dy, ch_dz);
			Tube driftTubeShape(rtubemin, rtubemax, chBox.z());
			Tube ShellTubeShape(rshellmin, rshellmax, chBox.z());
			Volume driftShellVolume(x_shell.nameStr(), ShellTubeShape, oddd.material(x_shell.materialStr()));
			driftShellVolume.setVisAttributes(oddd, x_shell.visStr());
			int ntubes = (chBox.x() - rshellmax) / rshellmax;
			double ygap = 2 * (ch_dy - 2 * (sqrt(3) * rshellmax + rshellmax));

			// Place the tubes and the chambers
			for (int j = 0; j < x_l.nz(); j++)
			{

				// place the chambers once in side A and once in side B
				chamberName = _toString((int)chamberNum, "chamber%d");
				double z = (j % 2) ? za : zb;
				Volume chVolume(chamberName, chBox, oddd.air());
				chVolume.setVisAttributes(oddd, x_ch.visStr());
				DetElement chamberElement(barrelMuonDetector, chamberName, chamberNum);
				PlacedVolume pv = barrelMuonVolume.placeVolume(chVolume, Transform3D(RotationZ(i * phistep), Position(x, y, z)));
				pv.addPhysVolID(x_ch.nameStr(), chamberNum++);
				chamberElement.setPlacement(pv);

				for (int nl = 0; nl < x_gas.repeat(); nl++)
				{

					// put the tubes (sensitive gas+shell) in two multilayers with ygap distance between
					double yt = -chBox.y() + rshellmax + nl * sqrt(3) * rshellmax;
					if (nl >= 0.5 * x_gas.repeat())
						yt += ygap;

					for (int t = 0; t < ntubes; t++)
					{

						tubeName = _toString((int)tubeNum, "tube%d");
						Volume driftTubeVolume(tubeName, driftTubeShape, oddd.material(x_gas.materialStr()));
						driftTubeVolume.setVisAttributes(oddd, x_gas.visStr());
						driftTubeVolume.setSensitiveDetector(sens);
						PlacedVolume pvtube = chVolume.placeVolume(driftTubeVolume, Transform3D(Position(-chBox.x() + (2 * t + 1) * rshellmax + (nl % 2) * rshellmax, yt, 0)));
						chVolume.placeVolume(driftShellVolume, Transform3D(Position(-chBox.x() + (2 * t + 1) * rshellmax + (nl % 2) * rshellmax, yt, 0)));
						pvtube.addPhysVolID("tube", tubeNum++);
						DetElement sensTubeEl(chamberElement, tubeName, tubeNum);
						sensTubeEl.setPlacement(pvtube);
					}
				}

				// move along z for the next chambers
				(j % 2) ? za -= zstep : zb += zstep;
			}
		}
	}

	// visualize the barrel cylinder
	barrelMuonVolume.setVisAttributes(oddd, x_det.visStr());

	// Place Volume
	Volume motherVolume = oddd.pickMotherVolume(barrelMuonDetector);
	Position translation(0., 0., x_det_dim.z());

	PlacedVolume placedMuonBarrel = motherVolume.placeVolume(barrelMuonVolume, translation);
	placedMuonBarrel.addPhysVolID("system", x_det.id());
	barrelMuonDetector.setPlacement(placedMuonBarrel);
	return barrelMuonDetector;
}

DECLARE_DETELEMENT(ODDMuonBarrel, create_element)
