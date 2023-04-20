#include "Yields.h"


Yields::Yields():
m_yield(0),
m_statistical(0),
m_systematic()
{}
Yields::Yields(double yield,double statistical):
m_yield(yield),
m_statistical(statistical)
{

m_systematic.push_back(0);
m_systematic.push_back(0);

}
Yields::~Yields()
{

  m_systematic.clear();

}
