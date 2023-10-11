#!/bin/bash

input=$1
output=$2

acts-install/bin/ActsAnalysisMaterialComposition \
      -i $input \
      -o $output/material_composition.root \
      --config ci/composition_config.json \
      -s

ci/make_detector_plot.py ci/composition_config.json $output/layout.pdf
