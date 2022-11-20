///////////////////////////////////////////////////////////
/// \file plot_macro.cpp
/// \brief ROOT macro to generate a plot from the
///        1D_path_integral.cpp output.
///
///  ROOT macro to generate a plot from the
///  1D_path_integral.cpp output of the method
///  Metropolis::ComputeEnergyEstimators. Various graphical
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
#include "TGraph.h"
#include "TGraphErrors.h"

int plot_macro()
{
#include "SETTINGS.h"

  // Initialize a canvas and a pad and set their attributes
  TCanvas* c = new TCanvas("c", "Data and Exact Solution", 800, 800);
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
  auto exact = new TGraph();

  // Add the exact line corresponding to constant DeltaE=1
  exact->AddPoint(-0.5, 1.);
  exact->AddPoint((double)N * 0.4 * a, 1.);
  exact->SetName("exact");

  // Open the data file
  ifstream infile(output_name);
  double time = 0., result = 0., error = 0., waste = 0.;

  // If open, extract all simulated data from the file
  if (infile) {
    infile.ignore(500, '\n');
    for (int i = 0; i < (int)((double)N * 0.4); i++) {
      infile >> time >> waste >> waste >> result >> error;
      montecarlo->AddPoint(time, result);
      montecarlo->SetPointError(i, 0., error);
    }
  } else {
    cout << "ERROR: \"output_file.dat\" not found\n";
    return 1;
  }

  // Tune some graphical aspects
  montecarlo->SetTitle("1D Anharmonic Oscillator, Simple Propagator");
  montecarlo->GetXaxis()->SetLimits(-0.2, (double)N * 0.34 * a);
  montecarlo->GetXaxis()->SetRangeUser(-0.2, (double)N * 0.34 * a);
  montecarlo->GetXaxis()->SetLabelSize(0.04);
  montecarlo->GetXaxis()->SetTitleSize(0.04);
  montecarlo->GetXaxis()->SetLabelOffset(0.008);
  montecarlo->GetXaxis()->SetTitleOffset(1.0);
  montecarlo->GetXaxis()->SetTitle("t");
  montecarlo->GetXaxis()->CenterTitle(true);
  montecarlo->GetYaxis()->SetTitle("#DeltaE(t)");
  montecarlo->GetYaxis()->SetLimits(0., 2.);
  montecarlo->GetYaxis()->SetRangeUser(0., 2.);
  montecarlo->GetYaxis()->SetLabelSize(0.04);
  montecarlo->GetYaxis()->SetTitleSize(0.04);
  montecarlo->GetYaxis()->SetLabelOffset(0.012);
  montecarlo->GetYaxis()->SetTitleOffset(2.);
  montecarlo->GetYaxis()->CenterTitle(true);
  montecarlo->SetName("montecarlo");
  montecarlo->SetMarkerColor(4);
  montecarlo->SetMarkerStyle(20);

  // Draw simulated data as points
  montecarlo->Draw("AP");

  // Draw the exact asymptotic function
  exact->SetLineColor(2);
  exact->Draw("lsame");

  // Initialize and draw a legend
  auto legend = new TLegend(0.15, 0.85, 0.45, 0.7);
  legend->AddEntry("exact", "Exact #DeltaE=1", "l");
  legend->AddEntry("montecarlo", "Montecarlo data", "ep");
  legend->SetTextSizePixels(18);
  legend->Draw("same");

  // Save plot into a file
  c->Print("plot_ROOT.png");
  return 0;
}
