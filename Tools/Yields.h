#ifndef YIELDS_H
#define YIELDS_H

#include <iostream>

class Yields {

 public:
 
  Yields();
  Yields(double yield,double statistical);
  ~Yields();
  
  inline double GetYield(){return m_yield;}
  inline double GetStatistical(){return m_statistical;}
  inline std::vector<double> GetSystematic(){return m_systematic;}
  inline void SetSystematic(std::vector<double> syst){m_systematic=syst;}

 private:
 
 double m_yield;
 double m_statistical;
 std::vector<double> m_systematic;


};

#endif // YIELDS_H
