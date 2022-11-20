///////////////////////////////////////////////////////////
/// \file plot_macro.cpp
/// \brief ROOT macro to generate a plot from the
///        VegasHarmOsc.cpp output.
///
///  ROOT macro to generate a plot from the
///  VegasHarmOsc.cpp output. This macro also performs
///  a fit to the function fitfunction, providing an
///  estimate for the ground energy E_0. Various graphical
///  options are implemented and tuned in the code. Type:\n
///  $ ROOT plot_macro.cpp\n
///  to run and generate the plot "plot_ROOT.png".
///
///////////////////////////////////////////////////////////
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"

///////////////////////////////////////////////////////////
/// Fit function
///
/// Function to fit data to: it represents the exact
/// asymptotic result from standard quantum mechanics
/// \param x pointer to the position x, where x[0]=x
/// \param par parameters of the function: par[0]=E_0, par[1]=T,
///        where E_0 is the ground state energy and T the propagation
///        time
///  \return \f[\left|\left\langle \psi |
///  E_{0}\right\rangle\right|^{2}e^{-E_{0}T}=\frac{1}{\sqrt{\pi}}e^{-x^{2}-E_{0}T}\f] \see
///  physical_params
///
///////////////////////////////////////////////////////////
double fitfunction(double* x, double* par)
{
  return exp(-x[0] * x[0] - par[0] * par[1]) / sqrt(3.14159);
}

int plot_macro()
{
#include "SETTINGS.h"

  // Define the root 1D function to be used in plot and fit
  TF1* f = new TF1("fit", fitfunction, -x_loop_bound, x_loop_bound, 2);

  // Set some useful properties, initialize parameters and fix the propagation time
  // to its tabulated value
  f->SetName("myfit");
  f->SetLineColor(2);
  f->SetParameter(0, 0.5);
  f->FixParameter(1, time_bound);
  f->SetParLimits(0, 0., 1.);
  f->SetParNames("E_0", "T");

  // Initialize a canvas and a pad and set their attributes
  TCanvas* c = new TCanvas("c", "Data and Fit", 800, 800);
  gStyle->SetOptFit(1011);
  gPad->SetGrid(1, 1);
  TPad* mypad = new TPad("mypad", "mypad", 0., 0., 1., 1.);
  mypad->Draw();
  mypad->cd();
  mypad->Range(0, 0, 100, 100);
  mypad->SetLeftMargin(0.15);
  mypad->SetRightMargin(0.15);
  mypad->SetTopMargin(0.15);
  mypad->SetBottomMargin(0.15);

  // Define the object containing the simulated data
  auto montecarlo = new TGraphErrors();

  // Open the data file
  ifstream infile("output_file.dat");
  string waste;
  double x = 0., result = 0., error = 0., exact = 0.;

  // If open, extract all simulated data from the file
  if (infile) {
    infile.ignore(500, '\n');
    for (int i = 0; i < (int)(2 * x_loop_bound / x_loop_step) + 1; i++) {
      infile >> x >> result >> error >> exact;
      montecarlo->AddPoint(x, result);
      montecarlo->SetPointError(i, 0., error);
    }
  } else {
    cout << "ERROR: \"output_file.dat\" not found\n";
    return 1;
  }

  // Tune some graphical aspects
  montecarlo->SetTitle("1D Harmonic Oscillator with Vegas");
  montecarlo->GetXaxis()->SetLimits(-x_loop_bound - 1., x_loop_bound + 1.);
  montecarlo->GetXaxis()->SetRangeUser(-x_loop_bound - 0.5, x_loop_bound + 0.5);
  montecarlo->GetXaxis()->SetLabelSize(0.04);
  montecarlo->GetXaxis()->SetTitleSize(0.04);
  montecarlo->GetXaxis()->SetLabelOffset(0.008);
  montecarlo->GetXaxis()->SetTitleOffset(1.0);
  montecarlo->GetXaxis()->SetTitle("x");
  montecarlo->GetXaxis()->CenterTitle(true);
  montecarlo->GetYaxis()->SetTitle("<x|exp(-HT)|x>");
  montecarlo->GetYaxis()->SetLimits(-1., 1.);
  montecarlo->GetYaxis()->SetRangeUser(0., 1.3 * (f->Eval(0.)));
  montecarlo->GetYaxis()->SetLabelSize(0.04);
  montecarlo->GetYaxis()->SetTitleSize(0.04);
  montecarlo->GetYaxis()->SetLabelOffset(0.012);
  montecarlo->GetYaxis()->SetTitleOffset(2.);
  montecarlo->GetYaxis()->CenterTitle(true);
  montecarlo->SetName("montecarlo");
  montecarlo->SetMarkerColor(4);
  montecarlo->SetMarkerStyle(20);

  // Perform the fit
  auto fitresult = montecarlo->Fit("myfit");
  // Draw the fit as a line and simulated data as points
  montecarlo->Draw("AP");

  // Draw the exact asymptotic function
  f->SetLineColor(3);
  f->SetParameter(0, 0.5);
  f->Draw("lsame");

  // Initialize and draw a legend
  auto legend = new TLegend(0.15, 0.85, 0.45, 0.7);
  legend->AddEntry("myfit", "#splitline{Exact asymptotic}{for E_0=1/2}", "l");
  legend->AddEntry("montecarlo", "Montecarlo data", "ep");
  legend->SetTextSizePixels(18);
  legend->Draw("same");

  // Save plot into a file
  c->Print("plot_ROOT.png");
  return 0;
}
