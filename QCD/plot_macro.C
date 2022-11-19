///////////////////////////////////////////////////////////
/// \file plot_macro.C
/// \brief ROOT macro to generate a plot from the
///        QCD_POST.C output.
///
///  ROOT macro to generate a plot from the
///  QCD_POST.C output of the method
///  Metropolis::ComputeRxTWilsonLoops which is stored in the file
///  RXT_potential_plot_file.dat. Various graphical
///  options are implemented and tuned in the code. Type:\n
///  $ root plot_macro.C\n
///  to run and generate the plot "QuarkPotential_root.png".
///
///////////////////////////////////////////////////////////
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"

/// Fit function
///
/// Function to fit data to: it represents a parametrization of the quark
/// potential, as a superposition of a short-distance coulomb profile and a
/// linear long-distance one, the latter representing confinement.
/// \param x pointer to the distance, where x[0]=r, the particles distance
/// \param par parameters of the function: par[0]=\f$\sigma\f$, par[1]=b, par[2]=c.
/// \return \f[V(r)=\sigma r-\frac{b}{r}+c\f]
double fitfunction(double* x, double* par)
{
  double value = -1.e8;
  if (x[0] != 0) value = par[0] * x[0] - par[1] / x[0] + par[2];
  return value;
}

/// Plot macro
///
/// Main function of the plot macro
int plot_macro()
{
  // Define the name of the input file containing data
  std::string filename = "RXT_potential_plot_file.dat";

  // Define the root 1D function to be used in plot and fit
  TF1* f = new TF1("fit", fitfunction, 0.1, 4.5, 3);

  // Set some useful properties, initialize parameters
  f->SetName("myfit");
  f->SetLineColor(2);
  f->SetParameters(5., 5., 0.);
  f->SetParLimits(0, -10., 10.);
  f->SetParLimits(1, -10., 10.);
  f->SetParLimits(2, -10., 10.);
  f->SetParNames("sigma", "b", "c");

  // Initialize a canvas and a pad and set their attributes
  TCanvas* c = new TCanvas("c", "Data and Fit", 800, 600);
  gStyle->SetOptFit(1011);

  // Define the object containing the simulated data
  auto gr = new TGraphErrors();

  // Open the data file
  std::ifstream infile(filename);
  std::string waste;
  int nR = 0;
  int n[4] = {0, 0, 0, 0};
  double r = 0., V = 0., V_error = 0.;

  // If open, extract all simulated data from the file
  if (infile) {
    infile.ignore(500, '\n');
    infile >> waste >> waste >> n[0] >> n[1] >> n[2] >> n[3];
    nR = (int)(min({n[0], n[1], n[2]}) / 2.) - 1;
    for (int i = 0; i < 14; i++) infile.ignore(500, '\n');
    for (int i = 0; i < nR; i++) {
      infile >> r >> V >> V_error;
      gr->AddPoint(r, V);
      gr->SetPointError(i, 0., V_error);
    }
  } else {
    std::cout << "ERROR: " << filename << " not found\n";
    return 1;
  }

  // Tune some graphical aspects
  gr->SetTitle("Quark Potential");
  gr->GetXaxis()->SetLimits(0., 5.);
  gr->GetXaxis()->SetRangeUser(0., 5.);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetLabelOffset(0.008);
  gr->GetXaxis()->SetTitleOffset(1.0);
  gr->GetXaxis()->SetTitle("r/a");
  gr->GetXaxis()->CenterTitle(true);
  gr->GetYaxis()->SetTitle("aV(r)");
  gr->GetYaxis()->SetLimits(-5., 15.);
  gr->GetYaxis()->SetRangeUser(0., 12.);
  gr->GetYaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetLabelOffset(0.012);
  gr->GetYaxis()->SetTitleOffset(1.0);
  gr->GetYaxis()->CenterTitle(true);
  gr->SetName("gr");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(20);

  // Perform the fit
  auto fitresult = gr->Fit("myfit", "S0"); // Don't draw the fit, yet
  // Draw simulated data as points
  gr->Draw("AP");

  // Set the fit function params to the fitted values and draw
  f->SetParameters(fitresult->Parameter(0), fitresult->Parameter(1), fitresult->Parameter(2));
  f->Draw("lsame");

  // Initialize and draw a legend
  auto legend = new TLegend(0.15, 0.85, 0.55, 0.65);
  legend->AddEntry("myfit", "Fit function V(r)", "l");
  legend->AddEntry("gr", "Montecarlo data", "ep");
  legend->SetTextSizePixels(25);
  legend->Draw("same");

  // Save plot into a file
  c->Print("QuarkPotential_root.png");
  return 0;
}
