#ifndef VARIABLEDISTR2D_H
#define VARIABLEDISTR2D_H

#include "VariableDistr.h"


class VariableDistr2D {

public:
  VariableDistr2D(VariableDistr *var_x,VariableDistr *var_y);

  ~VariableDistr2D();

  inline VariableDistr *GetVarX(){return m_var_x;}
  inline VariableDistr *GetVarY(){return m_var_y;}

private:
  VariableDistr *m_var_x;
  VariableDistr *m_var_y;

};


#endif // VARIABLEDISTR2D_H
