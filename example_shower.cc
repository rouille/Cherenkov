#include <iostream>
#include <vector>
#include <cmath>

#include "atmosphere.h"
#include "common.h"
#include "cherenkov.h"
#include "shower.h"
#include "conversion.h"

#include <TLegend.h>
#include <TRint.h>
#include <TRoot.h>
#include <TStyle.h>
#include <TGraphErrors.h>

static char gObjName[1024];
static int gObjNumber = 0;
static char * GetObjName() {sprintf(gObjName,"Shower %d",gObjNumber++); return gObjName;}


using namespace kPhysicalConstants;
using namespace std;



void Usage(string myName)
{
  cout << endl;
  cout << " Synopsis : " << endl;
  cout << myName << " <log(energy/[eV])> <zenith angle> <azimuth angle>" << endl << endl;

  cout << " Description :" << endl;
  cout << myName << " generates a shower produced by a cosmic ray of energy <log(energy/[eV])> and incoming direction <zenith angle>"
                 << " and <azimuth angle> and plots some of its characteristics" << endl;

  cout << endl;
  exit(0);
}



int main(int argc, char* argv[])
{
  // Command line
  if(argc != 4) Usage(argv[0]);
  double energy = pow(10,atof(argv[1]));
  double coord[2] = {atof(argv[2]),atof(argv[3])};

  // ROOT
  gROOT->SetStyle("Plain");
  gStyle->SetTitleFont(30,"TITLE");

  // Generate shower
  unsigned int NumberOfShower = 100000;
  unsigned int keep = 100;

  TCanvas * cFirstInteraction = new TCanvas(GetObjName(),"First Interaction",600,600);
  TH1F * hFirstInteraction = new TH1F(GetObjName(),"",50,0,40);
  hFirstInteraction->SetStats(0);

  vector< vector<double> > Ne_MC(keep);
  vector<double> T_MC;

  for(unsigned int i = 0; i < NumberOfShower; i++)
    {
      TShower * Shower = new TShower(energy,coord);
      double altitude = depth2altitude(Shower->GetT1()*X0);
      hFirstInteraction->Fill(altitude);

      if( i < keep )
        {
          vector<double> T_tmp, Ne_tmp;
          Shower->GenerateShower();
          Shower->GetLongitudinalProfile(T_tmp,Ne_tmp);
          Ne_MC[i].resize(T_tmp.size());
          for(unsigned j = 0; j < T_tmp.size(); j++) Ne_MC[i][j] = Ne_tmp[j];
          if( i == 0 ) T_MC = T_tmp;
        }

      delete Shower;
    }
  SetPlotAttributes(hFirstInteraction,"","Altitude of first interaction [km]");
  hFirstInteraction->Draw("");
  cFirstInteraction->Update();

  cFirstInteraction->SaveAs("X1.pdf");

  // Longitudinal profile
  TCanvas * cShower = new TCanvas(GetObjName(),"Shower",600,600);

  /* Individual showers */
  for(unsigned int i = 0; i < keep; i++)
    {
      TGraphErrors * gShower = new TGraphErrors(T_MC.size());
      for(unsigned int j = 0; j < T_MC.size(); j++) gShower->SetPoint(j,T_MC[j],log10(Ne_MC[i][j]));
      if( i == 0 )
        {
          SetPlotAttributes(gShower,0,40,"","t","N_{e}");
          gShower->GetYaxis()->SetRangeUser(0.,log10(energy)-8);
          gShower->Draw("AL");
        }
      else gShower->Draw("same L");
      cShower->Update();
    }

  /* Greisen */
  vector<double> T_Greisen = Bins(40, 0.1,40);
  vector<double> Ne_Greisen = Greisen(T_Greisen,energy);
  TGraphErrors * gGreisen = new TGraphErrors(T_Greisen.size());
  for(unsigned int i = 0; i < T_Greisen.size(); i++) gGreisen->SetPoint(i,T_Greisen[i],log10(Ne_Greisen[i]));
  gGreisen->SetMarkerStyle(kOpenCircle);
  gGreisen->SetMarkerSize(1.5);
  gGreisen->SetMarkerColor(kBlue);
  gGreisen->Draw("same P");
  cShower->Update();

  /* Mean */
  vector<double> Ne_Mean(T_MC.size());
  double mean = 0;
  for(unsigned int i = 0; i < T_MC.size(); i++)
    {
      mean = 0.;
      for(unsigned int j = 0; j < keep; j++) mean += Ne_MC[j][i];
      mean = mean*1./(1.*keep);
      Ne_Mean[i] = mean;
    }
  TGraphErrors * gMean = new TGraphErrors(T_MC.size());
  for(unsigned int i = 0; i < T_MC.size(); i++) gMean->SetPoint(i,T_MC[i],log10(Ne_Mean[i]));
  gMean->SetLineStyle(2);
  gMean->SetLineWidth(4);
  gMean->SetLineColor(kRed);
  gMean->Draw("same L");
  cShower->Update();

  /* Legend */
  char header[128]; sprintf(header,"Log(E_{0} / [eV]) = %d",(unsigned int)log10(energy));
  TLegend * legend = new TLegend(0.35,0.2,0.5,0.35);
  legend->SetFillColor(0); legend->SetBorderSize(0); legend->SetTextSize(0.035); legend->SetTextFont(132);
  legend->SetHeader(header);
  legend->AddEntry(gMean, "Mean", "L");
  legend->AddEntry(gGreisen, "Greisen", "P");
  legend->Draw();
  cShower->Update();


  cShower->SaveAs("shower.pdf");

  cout << "Program Finished Normally" << endl;
}
