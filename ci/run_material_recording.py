#!/usr/bin/env python3
import argparse
import acts, acts.examples
from pathlib import Path
import acts.examples.dd4hep
import acts.examples.geant4.dd4hep
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
import shutil
import tempfile
from acts.examples import (
    GaussianVertexGenerator,
    ParametricParticleGenerator,
    FixedMultiplicityGenerator,
    EventGenerator,
    RandomNumbers,
)
import os

u = acts.UnitConstants

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

parser = argparse.ArgumentParser(description="OpenDataDetector material recording")
parser.add_argument(
    "-n", "--events", type=int, default=100, help="Number of events to run"
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
    type=Path,
    default=Path.cwd(),
    help="Output directories. Default: $PWD",
)
args = parser.parse_args()
if args.jobs == -1: 
    args.jobs = multiprocessing.cpu_count()


def runMaterialRecording(seed: int, events: int, outputFile: Path):
    with tempfile.TemporaryDirectory() as outputDir:
        outputDir = Path(outputDir)
        oddDir = Path(__file__).parent.parent

        tracksPerEvent = 10

        detector, trackingGeometry, decorators = getOpenDataDetector(oddDir)

        detectorConstructionFactory = (
            acts.examples.geant4.dd4hep.DDG4DetectorConstructionFactory(detector)
        )

        #  dd4hepSvc = acts.examples.dd4hep.DD4hepGeometryService(
        #  xmlFileNames=[str(oddDir/ "xml/OpenDataDetector.xml")]
        #  )
        #  g4geo = acts.examples.geant4.dd4hep.DDG4DetectorConstruction(dd4hepSvc)

        s = acts.examples.Sequencer(events=events, numThreads=1)

        rnd = RandomNumbers(seed=seed)

        s = s or acts.examples.Sequencer(events=events, numThreads=1)

        evGen = EventGenerator(
            level=acts.logging.INFO,
            generators=[
                EventGenerator.Generator(
                    multiplicity=FixedMultiplicityGenerator(n=1),
                    vertex=GaussianVertexGenerator(
                        stddev=acts.Vector4(0, 0, 0, 0),
                        mean=acts.Vector4(0, 0, 0, 0),
                    ),
                    particles=ParametricParticleGenerator(
                        p=(1 * u.GeV, 10 * u.GeV),
                        eta=(-4, 4),
                        numParticles=tracksPerEvent,
                    ),
                )
            ],
            outputParticles="particles_initial",
            outputVertices="vertices_initial",
            randomNumbers=rnd,
        )

        s.addReader(evGen)

        g4Alg = acts.examples.geant4.Geant4MaterialRecording(
            level=acts.logging.INFO,
            detectorConstructionFactory=detectorConstructionFactory,
            randomNumbers=rnd,
            inputParticles=evGen.config.outputParticles,
            outputMaterialTracks="material_tracks",
        )

        #  g4AlgCfg = acts.examples.geant4.materialRecordingConfig(
        #  level=acts.logging.INFO,
        #  detector=g4geo,
        #  inputParticles=evGen.config.outputParticles,
        #  outputMaterialTracks="material_tracks",
        #  )

        #  g4AlgCfg.detectorConstruction = g4geo

        #  g4Alg = acts.examples.geant4.Geant4Simulation(
        #  level=acts.logging.INFO, config=g4AlgCfg
        #  )

        s.addAlgorithm(g4Alg)

        outfile = outputDir/ "geant4_material_tracks.root"
        s.addWriter(
            acts.examples.RootMaterialTrackWriter(
                prePostStep=True,
                recalculateTotals=True,
                inputMaterialTracks="material_tracks",
                filePath=str(outfile),
                level=acts.logging.INFO,
            )
        )

        s.run()

        del s

        shutil.copyfile(outfile, outputFile)


args.output.mkdir(exist_ok=True, parents=True)

with ProcessPoolExecutor(args.jobs) as ex:
    futures = []

    events_per_proc = args.events // args.jobs
    #  print(events_per_proc)
    events_per_proc = [events_per_proc] * args.jobs
    #  print(sum(events_per_proc))
    events_per_proc[-1] += args.events - sum(events_per_proc)
    #  print(sum(events_per_proc))

    for i, events in enumerate(events_per_proc):
        if events == 0:
            continue
        #  out = Path(f"mp_output/{i}")
        #  out.mkdir(exist_ok=True, parents=True)
        futures.append(
            ex.submit(
                runMaterialRecording,
                1247 * (i + 1),
                events,
                args.output / f"geant4_material_tracks_{i:>02d}.root",
            )
        )

    n = 0

    for f in as_completed(futures):
        n+=1
        print(n, " / ", len(futures))
        f.result()
