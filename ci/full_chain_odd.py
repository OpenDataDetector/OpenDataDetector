#!/usr/bin/env python3
import argparse
import pathlib, acts, acts.examples
import acts.examples.dd4hep

def getOpenDataDetector(odd_dir, mdecorator=None):
    import acts.examples.dd4hep
    dd4hepConfig = acts.examples.dd4hep.DD4hepGeometryService.Config(
        xmlFileNames=[str(odd_dir / "xml/OpenDataDetector.xml")]
    )
    detector = acts.examples.dd4hep.DD4hepDetector()

    config = acts.MaterialMapJsonConverter.Config()
    if mdecorator is None:
        mdecorator = acts.JsonMaterialDecorator(
            rConfig=config,
            jFileName=str(odd_dir / "config/odd-material-mapping-config.json"),
            level=acts.logging.WARNING,
        )

    trackingGeometry, deco = detector.finalize(dd4hepConfig, mdecorator)

    return detector, trackingGeometry, deco

parser = argparse.ArgumentParser(description="OpenDataDetector full chain example")
parser.add_argument(
    "-n", "--events", type=int, default=100, help="Number of events to run"
)
parser.add_argument(
    "-s", "--skip", type=int, default=0, help="Number of events to skip"
)
parser.add_argument(
    "-j",
    "--jobs",
    type=int,
    default=-1,
    help="Number of threads to use. Default: -1 i.e. number of cores",
)
parser.add_argument(
    "-o",
    "--output",
    type=pathlib.Path,
    default=pathlib.Path.cwd(),
    help="Output directories. Default: $PWD",
)
args = parser.parse_args()


u = acts.UnitConstants
outputDir = args.output

oddDir = pathlib.Path(__file__).parent.parent

oddMaterialMap = oddDir / "data/odd-material-maps.root"
oddDigiConfig = oddDir / "config/odd-digi-smearing-config.json"
oddSeedingSel = oddDir / "config/odd-seeding-config.json"
oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

detector, trackingGeometry, decorators = getOpenDataDetector(oddDir, mdecorator=oddMaterialDeco)
field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
rnd = acts.examples.RandomNumbers(seed=42)

from particle_gun import addParticleGun, MomentumConfig, EtaConfig, ParticleConfig
from fatras import addFatras
from digitization import addDigitization
from seeding import addSeeding, TruthSeedRanges
from ckf_tracks import addCKFTracks, CKFPerformanceConfig
from vertex_fitting import addVertexFitting, VertexFinder

s = acts.examples.Sequencer(events=args.events, numThreads=args.jobs, skip=args.skip)
s = addParticleGun(
    s,
    MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, True),
    EtaConfig(-3.0, 3.0, True),
    ParticleConfig(1, acts.PdgParticle.eMuon, True),
    rnd=rnd,
)
s = addFatras(
    s,
    trackingGeometry,
    field,
    outputDirRoot=outputDir,
    rnd=rnd,
)
s = addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=oddDigiConfig,
    outputDirRoot=outputDir,
    rnd=rnd,
)
s = addSeeding(
    s,
    trackingGeometry,
    field,
    TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-2.7, 2.7), nHits=(9, None)),
    geoSelectionConfigFile=oddSeedingSel,
    outputDirRoot=outputDir,
    initialVarInflation=[100, 100, 100, 100, 100, 100],
)
s = addCKFTracks(
    s,
    trackingGeometry,
    field,
    CKFPerformanceConfig(ptMin=400.0 * u.MeV, nMeasurementsMin=6),
    outputDirRoot=outputDir,
)
s = addVertexFitting(
    s,
    field,
    vertexFinder=VertexFinder.Truth,
    outputDirRoot=outputDir,
)

s.run()
