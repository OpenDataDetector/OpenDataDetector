// Open Data Dector project
//
// (c) 2021 CERN for the benefit of the ODD project
//
// Mozilla Public License Version 2.0

#pragma once

#include "DD4hep/DetFactoryHelper.h"

using namespace std;
using namespace dd4hep;

namespace ODDHelper {

template <typename T>
T& ensureExtension(DetElement& elt) {
  T* ext = elt.extension<T>(false);
  if (ext == nullptr) {
    ext = new T();
  }
  elt.addExtension<T>(ext);
  return *ext;
}

};  // namespace ODDHelper
