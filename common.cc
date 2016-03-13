#include "common.h"

#include <iostream>
#include <cmath>
#include <cstring>
#include <sys/stat.h>


bool CheckFile(string fileName)
{
  char fileNameCopy[100];
  struct stat fileStat;
  strcpy(fileNameCopy,fileName.c_str());
  int status = lstat(fileNameCopy,&fileStat);

  if( status == -1 )
  {   
    cout << "-----------------------------------------------------" << endl;
    cout << "File " << fileName << " does not exist." << endl;
    cout << "-----------------------------------------------------" << endl;
    return false;
  }   
  cout << "-----------------------------------------------------" << endl;
  cout << "Checking " << fileName << " (" << fileStat.st_size << " bytes)." << endl;
  if( fileStat.st_size == 0 ) 
  {   
    cout << "file size = 0 !!" << endl << endl;
    return false;
  }   
  cout << "-----------------------------------------------------" << endl;
  return true;
}



void SetPlotAttributes(TGraphErrors * Graph, double xMin, double xMax, string name, string Xaxis, string Yaxis)
{
  Graph->SetTitle(name.c_str());
  // X axis options
  Graph->GetHistogram()->SetXTitle(Xaxis.c_str());
  Graph->GetHistogram()->SetAxisRange(xMin,xMax);
  Graph->GetXaxis()->SetTitleFont(132);
  Graph->GetXaxis()->SetLabelSize(0.04);
  Graph->GetXaxis()->SetLabelFont(132);
  Graph->GetXaxis()->SetTitleSize(0.04);
  Graph->GetXaxis()->SetTitleOffset(1.1);
  // Y axis options
  Graph->GetHistogram()->SetYTitle(Yaxis.c_str());
  Graph->GetYaxis()->SetTitleFont(132);
  Graph->GetYaxis()->SetLabelSize(0.04);
  Graph->GetYaxis()->SetLabelFont(132);
  Graph->GetYaxis()->SetTitleSize(0.04);
  Graph->GetYaxis()->SetTitleOffset(1.3);
  // General
  Graph->SetMarkerStyle(kFullCircle);
  Graph->SetMarkerColor(kBlack);
  Graph->SetMarkerSize(1.);
}



void SetPlotAttributes(TH1F * Histo, string name, string Xaxis)
{
  Histo->SetTitle(name.c_str());
  // X axis options
  Histo->SetXTitle(Xaxis.c_str());
  Histo->GetXaxis()->SetTitleFont(132);
  Histo->GetXaxis()->SetLabelSize(0.03);
  Histo->GetXaxis()->SetLabelFont(132);
  Histo->GetXaxis()->SetTitleSize(0.035);
  Histo->GetXaxis()->SetTitleOffset(1.1);
  // Y axis options
  Histo->SetYTitle("Count");
  Histo->GetYaxis()->SetTitleFont(132);
  Histo->GetYaxis()->SetLabelSize(0.03);
  Histo->GetYaxis()->SetLabelFont(132);
  Histo->GetYaxis()->SetTitleSize(0.035);
  Histo->GetYaxis()->SetTitleOffset(1.4);
  // General
  Histo->SetStats(0);
}



vector<double> Bins(unsigned int size,double min,double max,bool logarithmic)
{
  vector<double> bins(size);
  vector<double> array(size);
  for(unsigned int i = 0; i < size; i++) array[i] = i/(size-1.);
  if( !logarithmic ) for(unsigned int i = 0; i < size; i++) bins[i] = min+(max-min)*array[i];
  else for(unsigned int i = 0; i < size; i++) bins[i] = pow(10,log10(min)+(log10(max)-log10(min))*array[i]);

  return bins;
}



double Integrate_nc5(const vector<double> & x, const vector<double> & y)
{
  // This is a five points Newton-Cotes (Bode's formula) integrator
  // See NumRec for details
  // We assume that the data is regularly gridded
  unsigned int npts = x.size();
  double h = x[1]-x[0];
  vector<unsigned int> ii; 
  unsigned int nbii = (unsigned int)floor((npts-1.)/4);
  unsigned int rest = (npts-1)-nbii*4;  
  unsigned int nbii2;
  if(rest == 1 || rest == 2) nbii2 = nbii-1;
  else nbii2 = nbii;
  ii.resize(nbii2);
  for(unsigned int i = 0; i < nbii2; i++) ii[i] = (i+1)*4;

  double integral = 0;
  for(unsigned int i = 0; i < nbii2; i++) integral += 2.*h*(7.*(y[ii[i]-4]+y[ii[i]])+32.*(y[ii[i]-3]+y[ii[i]-1])+12.*y[ii[i]-2])/45.;

  if( rest+1 == 2 ) 
  {
    // decoupage 4-3
    unsigned int shift = npts-2;
    // 4 points
    integral += 3*h*(y[shift-4]+3*y[shift-3]+3*y[shift-2]+y[shift-1])/8.;
    // 3 points
    integral += h*(y[npts-3]+4*y[npts-2]+y[npts-1])/3.;
  }
  else if( rest+1 == 3 )
  {
    // decoupage 4-4
    unsigned int shift = npts-3;
    // 4 points
    integral += 3*h*(y[shift-4]+3*y[shift-3]+3*y[shift-2]+y[shift-1])/8.;
    // 3 points
    shift = npts;
    integral += 3*h*(y[shift-4]+3*y[shift-3]+3*y[shift-2]+y[shift-1])/8.;
  }
  else if(rest+1 == 4)
  {
    integral += 3*h*(y[npts-4]+3*y[npts-3]+3*y[npts-2]+y[npts-1])/8.;
  }
  return integral;
}



vector<double> Interpol(const vector<double>& x, const vector<double>& y, const vector<double>& u)
{
  unsigned int size = u.size();
  unsigned int xsize = x.size();
  vector<double> v(size);
  unsigned int k, klow, khigh;
  for( unsigned int i=0;i<size;i++ )
  {   
    klow = 0;
    khigh = xsize-1;
    while(khigh-klow>1)
    {   
      k = (khigh+klow)/2.;
      if( x[k]>u[i] ) khigh = k;
      else klow = k;
    }   
    v[i] = y[klow]+((y[klow]-y[khigh])/(x[klow]-x[khigh]))*(u[i]-x[klow]);
  }   

  return v;
}


double Interpol(const vector<double>& x, const vector<double>& y, double u)
{
  vector<double> v(1);
  vector<double> uu(1);
  uu[0] = u;
  v = Interpol(x,y,uu);

  return v[0];
}

