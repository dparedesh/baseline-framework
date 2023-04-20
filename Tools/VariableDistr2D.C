#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>

using namespace std;

#include "VariableDistr.h"
#include "VariableDistr2D.h"

VariableDistr2D::VariableDistr2D(VariableDistr *var_x,VariableDistr *var_y):
  m_var_x(var_x),
  m_var_y(var_y)
{

}

VariableDistr2D::~VariableDistr2D(){


   

}
