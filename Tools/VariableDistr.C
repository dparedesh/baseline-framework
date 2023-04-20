#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>

using namespace std;


#include "VariableDistr.h"






VariableDistr::VariableDistr(const TString &variable,const TString &title,
                             const TString x_label, const TString y_label,
                             const Int_t nbBins,
                             const float xMin,const float xMax,bool yield,bool logX,bool logY) :
  m_variable(variable),
  m_title(title),
  m_nbBins(nbBins),
  m_xMin(xMin),
  m_xMax(xMax),
  m_yield(yield),
  m_x_label(x_label),
  m_y_label(y_label),
  m_logX(logX),
  m_logY(logY),
  m_array(false)
{
  

   float binning[]={};
   std::copy(binning,binning,m_binning);

 

}

VariableDistr::VariableDistr(const TString &variable,const TString &title,
                             const TString x_label, const TString y_label,
                             Int_t nbBins,
                             double binning[]) :
  m_variable(variable),
  m_title(title),
  m_nbBins(nbBins),
  m_yield(false),
  m_logX(false),
  m_logY(false),
  m_x_label(x_label),
  m_y_label(y_label)
{
  m_xMin=binning[0];
  m_xMax=binning[nbBins];

  int nEle=nbBins+1;

  m_binning = new double[nEle];
  
  memcpy(m_binning, binning, nEle*sizeof(double));

  m_array=true;
}


VariableDistr::~VariableDistr(){

}

TString VariableDistr::GetXLabel(){

  return m_x_label;
}
TString VariableDistr::GetYLabel(){

  float bin=(m_xMax-m_xMin)/m_nbBins;

  TString label=m_y_label;
  if (m_y_label.Contains("/")) label=m_y_label+" "+toString(bin);

  if (m_array){
        label=m_y_label+"Bin width";
        return label;
  }

  if (m_x_label.Contains("GeV")) label=label+" GeV";


return label;
}
string VariableDistr::toString(float number){
    ostringstream buff;
    buff<<number;
    return buff.str();
}
