// Open Data Dector project
//
// (c) 2021 CERN for the benefit of the ODD project
//
// Mozilla Public License Version 2.0

#pragma once

#include "DD4hep/DetFactoryHelper.h"

/// This method assemples a trapezoidal module for the
/// strip detectors
///
/// @param odd is the top level detector
/// @param sens is the top level sensitive detector container
/// @param x_module is the xml component describing the module
///
/// It excpects `module_component` xml childs
///
// @return a pair for a template module assembly and detector element
struct ODDModuleHelper
{

    static std::pair<dd4hep::Assembly, dd4hep::DetElement> assembleTrapezoidalModule(
        dd4hep::Detector &oddd, dd4hep::SensitiveDetector &sens, const xml_comp_t &x_module);

    /// This method assemples a rectangular module for the
    /// pixel and strip detectors
    ///
    /// @param odd is the top level detector
    /// @param sens is the top level sensitive detector container
    /// @param x_module is the xml component describing the module
    /// @param ylength[in,out] is the maximal length in y of all components
    ///        to be used for stave building
    ///
    /// It excpects `module_component` xml childs
    ///
    // @return a pair for a template module assembly and detector element
    static std::pair<dd4hep::Assembly, dd4hep::DetElement> assembleRectangularModule(
        dd4hep::Detector &oddd, dd4hep::SensitiveDetector &sens, const xml_comp_t &x_module,
        double &ylength);
};
