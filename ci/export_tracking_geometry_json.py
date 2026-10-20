#!/usr/bin/env python3
import argparse
from pathlib import Path

import acts
import acts.examples.dd4hep
from acts.json import TrackingGeometryJsonConverter
from acts.examples.odd import getOpenDataDetector

parser = argparse.ArgumentParser(
    description="Export the gen3 ODD tracking geometry to JSON"
)
parser.add_argument(
    "--odd-dir",
    type=Path,
    default=Path(__file__).parent.parent,
    help="ODD source directory. Default: this checkout",
)
parser.add_argument(
    "-o",
    "--output",
    type=Path,
    default=Path("odd-gen3.json"),
    help="Output file, format chosen by extension. Default: odd-gen3.json",
)
args = parser.parse_args()

detector = getOpenDataDetector(odd_dir=args.odd_dir, gen3=True)
trackingGeometry = detector.trackingGeometry()
gctx = acts.GeometryContext.dangerouslyDefaultConstruct()

converter = TrackingGeometryJsonConverter()
converter.toFile(gctx, trackingGeometry, str(args.output))

# Make sure the written file can be read back
converter.fromFile(gctx, str(args.output))
print(f"Wrote {args.output}")
