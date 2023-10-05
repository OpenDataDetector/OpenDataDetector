import ROOT
import argparse
import numpy
parser = argparse.ArgumentParser(description='Analyse calo shower data')
parser.add_argument('--infiles', '-i', required=True, type=str, nargs='+', help='ROOT files to analyse')
parser.add_argument('-o', '--outfile', type=str, default="showerAnalysis.root", help='output file')
parser.add_argument('--endcap', action='store_true', help='Perform analysis for endcap instead of barrel')
parser.add_argument('--hcal', action='store_true', help='Perform analysis for HCal instead of ECal')
args = parser.parse_args()

#__________________________________________________________
def draw_text(lines, coordinates = [0.1,0.8,0.5,0.9], colour = 36, border = 1):
   text = ROOT.TPaveText(coordinates[0],
                    coordinates[1],
                    coordinates[2],
                    coordinates[3],"brNDC")
   text.SetFillColorAlpha(0,0)
   text.SetBorderSize(border)
   text.SetTextFont(62)
   for line in lines:
      text.AddText("#color["+str(colour)+"]{"+line+"}")
      text.Draw()
      ROOT.SetOwnership(text,False)
   return text

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
def run(inputlist, outname):
    if ".root" not in outname:
        outname+=".root"
    # Read input parameters and create graphs
    g_resolution = ROOT.TGraphErrors()
    g_linearity = ROOT.TGraphErrors()
    for ifile, filename in enumerate(inputlist):
        file_in = ROOT.TFile(filename, "READ")
        for res in file_in.results:
            eMC = res.enMC
            e_mean = res.en_mean
            e_meanErr = res.en_meanErr
            e_resolution = res.en_resolution
            e_resolutionErr = res.en_resolutionErr
            g_resolution.SetPoint(ifile, eMC, e_resolution)
            g_resolution.SetPointError(ifile, 0, e_resolutionErr)
            g_linearity.SetPoint(ifile, eMC, e_mean / eMC)
            g_linearity.SetPointError(ifile, 0, e_meanErr / eMC)

    # Fit energy resolution
    f_fitResolution = ROOT.TF1("resolution", "sqrt([0]*[0] + pow([1]/sqrt(x),2))", g_resolution.GetXaxis().GetXmin(), g_resolution.GetXaxis().GetXmax())
    f_fitResolution.SetParName(0,"const")
    f_fitResolution.SetParName(1,"sqrt")
    result = g_resolution.Fit(f_fitResolution, 'S')
    formula = "#frac{#sigma_{E}}{E} = " + str(round(result.Get().Parameter(0)*100,2))+"% #oplus #frac{"+str(round(result.Get().Parameter(1)*100,2))+"%}{#sqrt{E}}"
    constString = "const term: "+str(round(result.Get().Parameter(0),4))+" #pm "+str(round(result.Get().Error(0),4))
    samplingString = "sampling term: "+str(round(result.Get().Parameter(1),4))+" #pm "+str(round(result.Get().Error(1),4))
    print(formula)
    print(constString)
    print(samplingString)

    # Fit linearity
    f_fitLinearity = ROOT.TF1("inearity", "pol0", g_linearity.GetXaxis().GetXmin(), g_linearity.GetXaxis().GetXmax())
    resultLin = g_linearity.Fit(f_fitLinearity, 'S')
    formulaLin = "#sum_{cells}dE/E_{MC} = " + str(round(resultLin.Get().Parameter(0)*100,2))+"% #pm "+str(round(resultLin.Get().Error(0),4))
    print(formulaLin)

    # Draw
    c_res = prepare_single_canvas("resolution","Energy resolution")
    prepare_graph(g_resolution, "resolution", ";E_{MC} (GeV);#sigma_{E_{sim}}/#LTE_{sim}#GT", 9, 21)
    g_resolution.Draw("ape")
    if not args.endcap:
        eta = 0
    else:
        eta = 2
    if not args.hcal:
        calo_type = "ECal resolution, photons"
    else:
        calo_type = "HCal resolution, pions"

    draw_text([calo_type+" at #||{#eta}="+str(eta)], [0.55,0.93,0.88,0.98], 1, 0).SetTextSize(0.05)
    draw_text([constString, samplingString], [0.55,0.68,0.88,0.76], 1, 0).SetTextSize(0.04)
    draw_text([formula], [0.55,0.8,0.88,0.9], 1, 0).SetTextSize(0.05)
    c_lin = prepare_single_canvas("linearity","Energy linearity")
    prepare_graph(g_linearity, "linearity", ";E_{MC} (GeV);#LTE_{sim}#GT/E_{MC}", 9, 21)
    g_linearity.Draw("ape")
    draw_text([calo_type + " at #||{#eta}="+str(eta)], [0.55,0.93,0.88,0.98], 1, 0).SetTextSize(0.05)
    draw_text([formulaLin], [0.55,0.8,0.88,0.9], 1, 0).SetTextSize(0.05)

    # Store
    c_res.SaveAs(f"resolution_{outname[:-5]}.pdf")
    c_lin.SaveAs(f"linearity_{outname[:-5]}.pdf")
    outfile = ROOT.TFile(outname, "RECREATE")
    outfile.cd()
    c_res.Write()
    c_lin.Write()
    g_resolution.Write("resolution")
    g_linearity.Write("linearity")
    const = numpy.zeros(1, dtype=float)
    sampl = numpy.zeros(1, dtype=float)
    response = numpy.zeros(1, dtype=float)
    constErr = numpy.zeros(1, dtype=float)
    samplErr = numpy.zeros(1, dtype=float)
    responseErr = numpy.zeros(1, dtype=float)
    tree = ROOT.TTree("params", "Fit parameters")
    tree.Branch("const", const, "const/D");
    tree.Branch("sampl", sampl, "sampl/D");
    tree.Branch("response", response, "response/D");
    tree.Branch("constErr", constErr, "constErr/D");
    tree.Branch("samplErr", samplErr, "samplErr/D");
    tree.Branch("responseErr", responseErr, "samplErr/D");
    const[0] = result.Get().Parameter(0)
    sampl[0] = result.Get().Parameter(1)
    response[0] = resultLin.Get().Parameter(0)
    constErr[0] = result.Get().Error(0)
    samplErr[0] = result.Get().Error(1)
    responseErr[0] = resultLin.Get().Error(0)
    tree.Fill()
    tree.Write()
    outfile.Close()

if __name__ == "__main__":

    run(args.infiles, args.outfile)
