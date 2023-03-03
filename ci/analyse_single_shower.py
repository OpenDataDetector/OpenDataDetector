import ROOT
import argparse
from math import sqrt
import numpy
parser = argparse.ArgumentParser(description='Analyse calo shower data')
parser.add_argument('--infile', '-i', required=True, type=str, nargs='+', help='EDM4hep file to analyse')
parser.add_argument('-o', '--outfile', type=str, default="showerAnalysis.root", help='output file')
parser.add_argument('-n', '--ncpus', type=int, default=2, help='Number of CPUs to use in analysis')
parser.add_argument('--endcap', action='store_true', help='Perform analysis for endcap instead of barrel')
args = parser.parse_args()

ROOT.gSystem.Load("libedm4hep")
#__________________________________________________________
def run(inputlist, outname, ncpu, endcapInsteadOfBarrel):
    outname = outname
    if ".root" not in outname:
        outname+=".root"
    if endcapInsteadOfBarrel:
        collname = "Endcap"
    else:
        collname = "Barrel"
    ROOT.ROOT.EnableImplicitMT(ncpu)
    df = ROOT.RDataFrame("events", inputlist)
    print ("Initialization done")
    # Create histograms with energy distributions
    h_ecal = df\
                 .Define("edep","ROOT::VecOps::RVec<float> result; for (auto&p: ECal"+collname+"Collection) {result.push_back(p.energy);} return result;")\
                 .Define("sumEdep","std::accumulate(edep.begin(),edep.end(),0.)")\
                 .Histo1D("sumEdep")
    h_mc = df\
               .Define("eMC","ROOT::VecOps::RVec<float> result; for(auto& m:MCParticles){result.push_back(sqrt(m.momentum.x*m.momentum.x+m.momentum.y*m.momentum.y+m.momentum.z*m.momentum.z));} return result;")\
               .Define("gunMC","eMC[0]")\
               .Histo1D("gunMC")
    h_ratio = df\
                  .Define("edep","ROOT::VecOps::RVec<float> result; for (auto&p: ECal"+collname+"Collection) {result.push_back(p.energy);} return result;")\
                  .Define("sumEdep","std::accumulate(edep.begin(),edep.end(),0.)")\
                  .Define("eMC","ROOT::VecOps::RVec<float> result; for(auto& m:MCParticles){result.push_back(sqrt(m.momentum.x*m.momentum.x+m.momentum.y*m.momentum.y+m.momentum.z*m.momentum.z));} return result;")\
                  .Define("gunMC","eMC[0]")\
                  .Define("eratio","sumEdep/gunMC")\
                  .Histo1D("eratio")
    print(f"E distribution EMCal: <E>= {h_ecal.GetMean()}\t RMS= {h_ecal.GetRMS()}")
    print(f"E distribution MC particles: E_MC= {h_mc.GetMean()}\t RMS= {h_mc.GetRMS()}")
    print(f"sampling fraction calculated as <E>/E_MC: {h_ratio.GetMean()}")
    gun_mean = h_mc.GetMean()
    gun_rms = h_mc.GetRMS()

    # Fitting Gaussian part
    f_prefit = ROOT.TF1("firstGaus","gaus", h_ecal.GetMean() - 2. * h_ecal.GetRMS(), h_ecal.GetMean() + 2. * h_ecal.GetRMS())
    result_pre = h_ecal.Fit(f_prefit, "SRQN")
    f_fit = ROOT.TF1("finalGaus", "gaus", result_pre.Get().Parameter(1) - 2. * result_pre.Get().Parameter(2), result_pre.Get().Parameter(1) + 2. * result_pre.Get().Parameter(2) )
    result = h_ecal.Fit(f_fit, "SRQ")
    result_mean = result.Get().Parameter(1)
    result_meanError = result.Get().Error(1)
    result_resolution = result.Get().Parameter(2) / result.Get().Parameter(1)
    tmp_resolutionErrorSigma = result.Get().Error(2) / result.Get().Parameter(1)
    tmp_resolutionErrorMean = result.Get().Error(1) * result.Get().Parameter(2) / ( result.Get().Parameter(1) ** 2)
    result_resolutionError = sqrt( tmp_resolutionErrorSigma ** 2 +  tmp_resolutionErrorMean ** 2 )
    print("Fitting Gaussian to EMCal energy distribution...")
    print(f"\tmean energy <E>= {result_mean}\t +- {result_meanError}")
    print(f"\tresolution sigma(E)= {result_resolution}\t +- {result_resolutionError}")
    print(f"\tsampling fraction calculated as <E>/E_MC: {result_mean/gun_mean}\t +- {result_meanError/gun_mean}")

    # Store
    outfile = ROOT.TFile(outname, "RECREATE")
    outfile.cd()
    h_ecal.Write("energy_ecal")
    h_mc.Write("energy_MC0")
    h_ratio.Write("energy_ratio")
    store_mean = numpy.zeros(1, dtype=float)
    store_meanErr = numpy.zeros(1, dtype=float)
    store_resolution = numpy.zeros(1, dtype=float)
    store_resolutionErr = numpy.zeros(1, dtype=float)
    store_EMC = numpy.zeros(1, dtype=float)
    tree = ROOT.TTree("results", "Fit parameters and E_MC")
    tree.Branch("en_mean", store_mean, "store_mean/D");
    tree.Branch("en_meanErr", store_meanErr, "store_meanErr/D");
    tree.Branch("en_resolution", store_resolution, "store_resolution/D");
    tree.Branch("en_resolutionErr", store_resolutionErr, "store_resolutionErr/D");
    tree.Branch("enMC", store_EMC, "store_EMC/D");
    store_mean[0] = result_mean
    store_meanErr[0] = result_meanError
    store_resolution[0] = result_resolution
    store_resolutionErr[0] = result_resolutionError
    store_EMC[0] = gun_mean
    tree.Fill()
    tree.Write()
    outfile.Close()

if __name__ == "__main__":

    run(args.infile, args.outfile, args.ncpus, args.endcap)
