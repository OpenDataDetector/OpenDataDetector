#!/bin/bash

input=$1
output=$2

exec acts-install/bin/ActsAnalysisMaterialComposition \
      -i $input \
      -o $output \
      --config composition_config.json \
      -s
