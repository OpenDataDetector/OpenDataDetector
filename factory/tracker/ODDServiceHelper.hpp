// Open Data Dector project
//
// (c) 2021 CERN for the benefit of the ODD project
//
// Mozilla Public License Version 2.0

#pragma once

#include "DD4hep/DetFactoryHelper.h"

/// Helper method to build service routing for the Barrel type
///
/// @tparam volume_t the Type of volume (Volume, Assembly)
///
/// @param odd the top level detector
/// @param barrelVolume the volume to put the routing on
/// @param x_routing the xml description of the routing
/// @param layerR the layer radii to connect
template <typename volume_t>
void buildBarrelRouting(dd4hep::Detector& oddd, volume_t& barrelVolume,
                        const xml_comp_t& x_routing,
                        const std::vector<double>& layerR) {

  using namespace dd4hep;

  // Grab the cables & route them outwards
  unsigned int nphi = x_routing.nphi();

  double phiStep = 2 * M_PI / nphi;
  double phi0 = x_routing.phi0();
  double rmin = x_routing.rmin();
  double rmax = x_routing.rmax();
  double n = x_routing.number();

  // Loop over the layer routings
  for (unsigned int ib = 1; ib < layerR.size(); ++ib) {
    double gap = x_routing.gap();
    double clength = layerR[ib] - layerR[ib - 1] - 2. * gap;
    double rpos = 0.5 * (layerR[ib] + layerR[ib - 1]);

    Tube cable(rmin, rmax, 0.5 * clength);
    Volume cableVolume("Cable", cable,
                       oddd.material(x_routing.materialStr()));
    cableVolume.setVisAttributes(oddd, x_routing.visStr());

    for (int side = -1; side < 2; side += 2) {

      Assembly cableboxAssembly("CableBox");
      if (x_routing.hasChild(_U(box))) {
        // The box plate for the cables
        xml_comp_t x_box = x_routing.child(_U(box));
        Box box(x_box.dz(), n * ib * rmax, 0.5 * clength);
        Volume boxVolume("CableBand", box,
                         oddd.material(x_routing.materialStr()));
        boxVolume.setVisAttributes(oddd, x_box.visStr());

        cableboxAssembly.placeVolume(
            boxVolume, Position(side * (rmax + x_box.dz()), 0., 0.));
      }

      for (unsigned int icable = 0; icable < n * ib; ++icable) {
        // Place the pipe in the stave
        cableboxAssembly.placeVolume(
            cableVolume, Position(0., (-n * ib + 1 + 2 * icable) * rmax, 0.));
      }

      for (unsigned int iphi = 0; iphi < nphi; ++iphi) {
        // Calculate the phi
        double phi = phi0 + iphi * phiStep;

        // The layer position
        double xpos = rpos * cos(phi);
        double ypos = rpos * sin(phi);
        double zpos = side * x_routing.z_offset();

        // Place the pipe in the stave
        barrelVolume.placeVolume(
            cableboxAssembly,
            Transform3D(RotationZ(phi) * RotationY(0.5 * M_PI),
                        Position(xpos, ypos, zpos)));
      }
    }
  }
}

/// Helper method to build service routing for the Endcap type
///
/// @tparam volume_t the Type of volume (Volume, Assembly)
///
/// @param odd the top level detector
/// @param endcapVolume the volume to put the routing on
/// @param x_routing the xml description of the routing
/// @param endcapZ the layer z positions to connect
template <typename volume_t>
void buildEndcapRouting(dd4hep::Detector& oddd, volume_t& endcapVolume,
                        const xml_comp_t& x_routing,
                        const std::vector<double>& endcapZ) {
  
  using namespace dd4hep;
  
  // Grab the cables & route them outwards
  unsigned int nphi = x_routing.nphi();

  double phiStep = 2 * M_PI / nphi;
  double phi0 = x_routing.phi0();
  double rmin = x_routing.rmin();
  double rmax = x_routing.rmax();
  double r = x_routing.r();
  double n = x_routing.number();

  // Loop over the layer routings
  for (unsigned int iec = 1; iec < endcapZ.size(); ++iec) {
    double gap = x_routing.gap();
    double clength = std::abs(endcapZ[iec] - endcapZ[iec - 1]) - 2. * gap;

    Assembly cableboxAssembly("CableBox");
    if (x_routing.hasChild(_U(box))) {
      // The box plate for the cables
      xml_comp_t x_box = x_routing.child(_U(box));
      Box box(x_box.dz(), n * iec * rmax, 0.5 * clength);
      Volume boxVolume("CableBand", box,
                       oddd.material(x_routing.materialStr()));
      boxVolume.setVisAttributes(oddd, x_box.visStr());

      cableboxAssembly.placeVolume(boxVolume,
                                   Position(rmax + x_box.dz(), 0., 0.));
    }

    Tube cable(rmin, rmax, 0.5 * clength);
    Volume cableVolume("Cable", cable,
                       oddd.material(x_routing.materialStr()));
    cableVolume.setVisAttributes(oddd, x_routing.visStr());

    for (unsigned int icable = 0; icable < n * iec; ++icable) {
      // Place the pipe in the stave
      cableboxAssembly.placeVolume(
          cableVolume, Position(0., (-n * iec + 1 + 2 * icable) * rmax, 0.));
    }

    for (unsigned int iphi = 0; iphi < nphi; ++iphi) {
      // Calculate the phi
      double phi = phi0 + iphi * phiStep;

      // The layer position
      double xpos = r * cos(phi);
      double ypos = r * sin(phi);
      double zpos = 0.5 * (endcapZ[iec] + endcapZ[iec - 1]);

      // Place the pipe in the stave
      endcapVolume.placeVolume(
          cableboxAssembly,
          Transform3D(RotationZ(+phi), Position(xpos, ypos, zpos)));
    }
  }
}

/// Helper method to build a cyoindrical like passive structure
///
/// @tparam volume_t the Type of volume (Volume, Assembly)
///
/// @param odd the top level detector
/// @param endcapVolume the volume to put the routing on
/// @param x_mother_comp the xml description of teh mother component
/// @param layerR the layer radii contaienr to add the new one
template <typename volume_t>
void buildSupportCylinder(dd4hep::Detector& oddd, volume_t& motherVolume,
                          const xml_comp_t& x_mother_comp,
                          std::vector<double>& layerR) {
  
  using namespace dd4hep;
  
  size_t supportNum = 0;
  for (xml_coll_t sup(x_mother_comp, _U(support)); sup; ++sup, ++supportNum) {
    xml_comp_t x_support = sup;
    // Create the volume of the support structure
    std::string supportName = _toString((int)supportNum, "SupportCylinder%d");

    // Remember the layer radius if it is needed for operation
    if (x_support.hasChild(_Unicode(connector))) {
      layerR.push_back(0.5 * (x_support.rmin() + x_support.rmax()));
    }
    // If nz is not set to 0, build 2 symmetric ones
    for (int side = -1; side < x_support.nsides(); side += 2) {
      // Create the support volume
      Volume supportVolume(
          supportName, Tube(x_support.rmin(), x_support.rmax(), x_support.dz()),
          oddd.material(x_support.materialStr()));
      supportVolume.setVisAttributes(oddd, x_support.visStr());
      // Place the support structure
      motherVolume.placeVolume(supportVolume,
                               Position(0., 0., side * x_support.z_offset()));
    }
  }
}

/// Helper method to build a cyoindrical like passive structure
///
/// @tparam volume_t the Type of volume (Volume, Assembly)
///
/// @param odd the top level detector
/// @param endcapVolume the volume to put the routing on
/// @param x_mother_comp the xml description of teh mother component
/// @param layerR the layer radii contaienr to add the new one
template <typename volume_t>
void buildCoolingRings(dd4hep::Detector& oddd, volume_t& motherVolume,
                       const xml_comp_t& x_mother_comp) {
  
  using namespace dd4hep;
  
  size_t cringNum = 0;
  for (xml_coll_t cring(x_mother_comp, _Unicode(cooling_ring)); cring;
       ++cring, ++cringNum) {
    xml_comp_t x_cooling_ring = cring;

    double r = x_cooling_ring.r();
    double nPhi = x_cooling_ring.nphi();
    double phiStep = 2. * M_PI / nPhi;
    double zpos = x_cooling_ring.z_offset();
    double dz = 2 * (r * M_PI / nPhi - x_cooling_ring.gap());

    Volume coolingSegement(
        "CoolingRingSegment",
        Tube(x_cooling_ring.rmin(), x_cooling_ring.rmax(), dz),
        oddd.material(x_cooling_ring.materialStr()));
    coolingSegement.setVisAttributes(oddd, x_cooling_ring.visStr());

    // Create the segments around the ring
    for (unsigned int iphi = 0; iphi < nPhi; ++iphi) {
      // position & orientation
      double phi = iphi * phiStep;
      Position segementPos(r * cos(phi), r * sin(phi), zpos);
      // Place the support structure
      motherVolume.placeVolume(
          coolingSegement,
          Transform3D(RotationY(0.5 * M_PI) * RotationX(0.5 * M_PI - phi),
                      segementPos));
    }
  }
}
