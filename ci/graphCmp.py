import ROOT
import argparse
import numpy
parser = argparse.ArgumentParser(description='Analyse calo shower data')
parser.add_argument('--infiles', '-i', required=True, type=str, nargs='+', help='ROOT files to draw')
parser.add_argument('--legend', '-l', required=True, type=str, nargs='+', help='Labels in legend corresponding to input')
parser.add_argument('-o', '--outfile', type=str, default="comparePerformance.root", help='output file')
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
    for ifile, filename in enumerate(inputlist):
        file_in = ROOT.TFile(filename, "READ")
        gR = file_in.Get("resolution")
        gL = file_in.Get("linearity")
        ROOT.SetOwnership(gR, False)
        ROOT.SetOwnership(gL, False)
        g_resolutions.append(gR)
        g_linearities.append(gL)

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

if __name__ == "__main__":

    run(args.infiles, args.legend, args.outfile)
