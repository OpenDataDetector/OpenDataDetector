import ROOT
import argparse
import numpy
import sys
parser = argparse.ArgumentParser(description='Analyse calo shower data')
parser.add_argument('--infiles', '-i', required=True, type=str, nargs='+', help='ROOT files to draw')
parser.add_argument('--legend', '-l', required=True, type=str, nargs='+', help='Labels in legend corresponding to input')
parser.add_argument('-o', '--outfile', type=str, default="compareToReference.root", help='output file')
args = parser.parse_args()

#__________________________________________________________
def prepare_single_canvas(name, title):
   c = ROOT.TCanvas(name, title, 1200, 900)
   c.SetTopMargin(0.01)
   c.SetRightMargin(0.03)
   c.SetLeftMargin(0.15)
   c.SetBottomMargin(0.15)
   c.SetGridx()
   c.SetGridy()
   ROOT.SetOwnership(c,False)
   return c
#__________________________________________________________
def prepare_graph(graph, name, title, colour = 9, markerStyle = 21, factor = 1):
   # graph settings
   graph.SetTitle(title)
   graph.SetName(name)
   graph.SetMarkerStyle(markerStyle)
   graph.SetMarkerSize(1.4)
   graph.SetMarkerColor(colour)
   graph.SetLineColor(colour)
   for fnc in graph.GetListOfFunctions():
       fnc.SetLineColor(colour)
   # set Y axis
   graph.GetYaxis().SetTitleSize(0.06)
   graph.GetYaxis().SetTitleOffset(1.1)
   graph.GetYaxis().SetLabelSize(0.045)
   graph.GetYaxis().SetNdivisions(504)
   # set X axis
   graph.GetXaxis().SetTitleSize(0.07)
   graph.GetXaxis().SetTitleOffset(1.)
   graph.GetXaxis().SetLabelSize(0.05)
   graph.GetYaxis().SetNdivisions(506)



#__________________________________________________________
def run(inputlist, legendNames, outname):
    if ".root" not in outname:
        outname+=".root"
    # Read input parameters and create graphs
    g_resolutions = []
    g_linearities = []
    v_resSampl = []
    v_resSamplErr = []
    v_resConst = []
    v_resConstErr = []
    v_resResponse = []
    v_resResponseErr = []
    for ifile, filename in enumerate(inputlist):
        file_in = ROOT.TFile(filename, "READ")
        gR = file_in.Get("resolution")
        gL = file_in.Get("linearity")
        params = file_in.Get("params")
        ROOT.SetOwnership(gR, False)
        ROOT.SetOwnership(gL, False)
        g_resolutions.append(gR)
        g_linearities.append(gL)
        for p in params:
           v_resSampl.append(p.sampl)
           v_resSamplErr.append(p.samplErr)
           v_resConst.append(p.const)
           v_resConstErr.append(p.constErr)
           v_resResponse.append(p.response)
           v_resResponseErr.append(p.responseErr)

    # Draw
    listColours = [9, 1, 2, 3, 94, 6]
    listMarkers = [21, 20, 24, 25, 33, 27]
    c_res = prepare_single_canvas("resolution","Energy resolution")
    legendRes = ROOT.TLegend(0.7,0.7,0.9,0.9)
    for igraph, graph in enumerate(g_resolutions):
        prepare_graph(graph, "resolution", ";E_{MC} (GeV);#sigma_{E_{sim}}/#LTE_{sim}#GT", listColours[igraph], listMarkers[igraph])
        legendRes.AddEntry(graph,legendNames[igraph],"lep");
        if igraph ==0:
            graph.Draw("ape")
        else:
            graph.Draw("samep")
    legendRes.Draw()
    c_lin = prepare_single_canvas("linearity","Energy linearity")
    legendLin = ROOT.TLegend(0.7,0.7,0.9,0.9)
    for igraph, graph in enumerate(g_linearities):
        prepare_graph(graph, "linearity", ";E_{MC} (GeV);#LTE_{sim}#GT/E_{MC}", listColours[igraph], listMarkers[igraph])
        legendLin.AddEntry(graph,legendNames[igraph],"lep");
        if igraph ==0:
            graph.Draw("ape")
        else:
            graph.Draw("samep")
    legendLin.Draw()

    # Store
    c_res.SaveAs(f"resolution_{outname[:-5]}.pdf")
    c_lin.SaveAs(f"linearity_{outname[:-5]}.pdf")
    outfile = ROOT.TFile(outname, "RECREATE")
    outfile.cd()
    c_res.Write()
    c_lin.Write()
    outfile.Close()

    # Compare resolution (for 2 input files only, second one is used as reference)
    if len(inputlist)==2:
       if v_resSampl[0] + 2 * v_resSamplErr[0] < v_resSampl[1] -  2 * v_resSamplErr[1] or v_resSampl[0] - 2 * v_resSamplErr[0] > v_resSampl[1] + 2 * v_resSamplErr[1]:
          print(f'\033[91m\nResolution has changed.\nSampling term of new sample is {v_resSampl[0]}+-{v_resSamplErr[0]} vs reference {v_resSampl[1]}+-{v_resSamplErr[1]}.\nConstant term will not be checked.\n\033[00m')
          sys.exit(1)
       else:
          print(f'\033[32m\nSampling term of new sample is {v_resSampl[0]}+-{v_resSamplErr[0]} which is within the error bars wrt reference sample {v_resSampl[1]}+-{v_resSamplErr[1]}\n\033[00m')
       if v_resConst[0] + 2 * v_resConstErr[0] < v_resConst[1] - 2 * v_resConstErr[1] or v_resConst[0] - 2 * v_resConstErr[0] > v_resConst[1] + 2 * v_resConstErr[1]:
          print(f'\033[91m\nResolution has changed.\nConstant term of new sample is {v_resConst[0]}+-{v_resConstErr[0]} vs reference {v_resConst[1]}+-{v_resConstErr[1]}.\n\033[00m')
          sys.exit(1)
       else:
          print(f'\033[32m\nConstant term of new sample is {v_resConst[0]}+-{v_resConstErr[0]} which is within the error bars wrt reference sample {v_resConst[1]}+-{v_resConstErr[1]}.\n\033[00m')
       if v_resResponse[0] + 2 * v_resResponseErr[0] < v_resResponse[1] - 2 * v_resResponseErr[1] or v_resResponse[0] - 2 * v_resResponseErr[0] > v_resResponse[1] + 2 * v_resResponseErr[1]:
          print(f'\033[91m\nReponse of the detector has changed.\nResponse of new sample is {v_resResponse[0]}+-{v_resResponseErr[0]} vs reference {v_resResponse[1]}+-{v_resResponseErr[1]}.\n\033[00m')
          sys.exit(1)
       else:
          print(f'\033[32m\nResponse of the detector for new sample is {v_resResponse[0]}+-{v_resResponseErr[0]} which is within the error bars wrt reference sample {v_resResponse[1]}+-{v_resResponseErr[1]} within the 2 sigma.\n\033[00m')

if __name__ == "__main__":

    run(args.infiles, args.legend, args.outfile)
