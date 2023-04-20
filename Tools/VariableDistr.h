#ifndef VARIABLEDISTR_H
#define VARIABLEDISTR_H


#include "TString.h"

class VariableDistr {

public:
  //VariableDistr();
  //
  //

  VariableDistr(const TString &variable,const TString &title,const TString x_label, const TString y_label,
                const Int_t nbBins,const float xMin,const float xMax,bool yield=false,bool logX=false,bool logY=false);

  VariableDistr(const TString &variable,const TString &title, const TString x_label, const TString y_label,
                const Int_t nbBins,double binning[]);

  ~VariableDistr();

  inline TString GetName(){return m_variable;}
  inline TString GetTitle(){return m_title;}
  inline int GetNbins(){return m_nbBins;}
  inline float GetLowAxis(){return m_xMin;}
  inline float GetUpAxis(){return m_xMax;}
  inline bool GetYield(){return m_yield;}
         TString GetXLabel();
         TString GetYLabel();
  string toString(float number);

  inline double* GetBinning(){return m_binning;};
  inline bool IsArray(){return m_array;};
  inline bool plotLogX(){return m_logX;};
  inline bool plotLogY(){return m_logY;};

private:
  TString m_variable;
  TString m_title;
  TString m_x_label;
  TString m_y_label;
  int m_nbBins;
  float m_xMin,m_xMax;
  bool m_yield;
  double *m_binning;
  bool m_array;
  bool m_logX;
  bool m_logY;

};


#endif // VARIABLEDISTR_H

