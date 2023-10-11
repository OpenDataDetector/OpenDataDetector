#!/bin/bash

input=$1
output=$2

acts-install/bin/ActsAnalysisMaterialComposition \
      -i $input \
      -o $output \
      --config ci/composition_config.json \
      -s

ci/make_material_plots.py ci/composition_config.json $output/layout.pdf
