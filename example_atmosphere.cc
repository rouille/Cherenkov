#include <iostream>
#include <vector>

#include "atmosphere.h"
#include "common.h"
#include "cherenkov.h"
#include "shower.h"

#include <TRint.h>
#include <TRoot.h>
#include <TStyle.h>

static char gObjName[1024];
static int gObjNumber = 0;
static char * GetObjName() {sprintf(gObjName,"Atmosphere %d",gObjNumber++); return gObjName;}



using namespace std;



void Usage(string myName)
{
  cout << endl;
  cout << " Synopsis : " << endl;
  cout << myName << " <atmospheric file>" << endl << endl;

  cout << " Description :" << endl;
  cout << myName << " plots the atmospheric parameters enclosed in <atmospheric file> and the Cherenkov condition in air" << endl;

  cout << endl;
  exit(0);
}



int main(int argc, char* argv[])
{
  // Command line
  if(argc != 2) Usage(argv[0]);
  string fileName = argv[1];
  if( !CheckFile(fileName) ) {cerr << "Exiting" << endl; exit(0);}

  // ROOT
  gROOT->SetStyle("Plain");
  gStyle->SetTitleFont(30,"TITLE");

  // Get the atmosphere
  vector<TAtmosphere> Atmosphere = GetAtmosphere(fileName);

  DECLARE_POINTER(double,altitude,Atmosphere,fAltitude);
  DECLARE_POINTER(double,density,Atmosphere,fDensity);
  DECLARE_POINTER(double,depth,Atmosphere,fDepth);
  DECLARE_POINTER(double,delta,Atmosphere,fDelta);

  // Plot
  TCanvas * cAtmosphere = new TCanvas(GetObjName(),"Atmoshere",2500,2500);
  cAtmosphere->Divide(2,2);

  /* altitude vs density */
  cAtmosphere->cd(1);
  TGraphErrors * gDensity = new TGraphErrors(Atmosphere.size(),density,altitude);
  SetPlotAttributes(gDensity,0.,2.,"","#rho [g . cm^{-3}]","Altitude [km]");
  gDensity->SetLineWidth(2);
  gDensity->Draw("AL");
  cAtmosphere->Update();

  /* altitude vs depth */
  cAtmosphere->cd(2);
  TGraphErrors * gDepth = new TGraphErrors(Atmosphere.size(),depth,altitude);
  SetPlotAttributes(gDepth,0.,1200.,"","Atmospheric depth [g . cm^{-2}]","Altitude [km]");
  gDepth->SetLineWidth(2);
  gDepth->Draw("AL");
  cAtmosphere->Update();

  /* altitude vs refractive index */
  cAtmosphere->cd(3);
  TGraphErrors * gDelta = new TGraphErrors(Atmosphere.size(),delta,altitude);
  SetPlotAttributes(gDelta,0.,0.0003,"","#delta","Altitude [km]");
  gDelta->SetLineWidth(2);
  gDelta->Draw("AL");
  cAtmosphere->Update();

  /* altitude vs Cherenkov energy threshold */
  cAtmosphere->cd(4);
  double * Eth = CherenkovEnergyThreshold(Atmosphere.size(),delta);
  cAtmosphere->GetPad(4)->SetLogx();
  TGraphErrors * gCherenkov = new TGraphErrors(Atmosphere.size(),Eth,altitude);
  SetPlotAttributes(gCherenkov,1,150000,"","Cherenkov Energy Threshold [MeV]","Altitude [km]");
  gCherenkov->Draw("AP");
  cAtmosphere->Update();

  cAtmosphere->SaveAs("atmosphere.pdf");

  delete altitude;
  delete density;
  delete depth;
  delete delta;
  delete Eth;

  cout << "Program Finished Normally" << endl;
}
