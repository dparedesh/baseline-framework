#include "EventCut.h"

EventCut::EventCut(TString cut,TString name):
// m_cut("("+cut+")"),
 m_cut(""),
 m_name(name),
 m_min(-99),
 m_max(-99)
{

  if (cut!="") m_cut="("+cut+")";


}

EventCut::~EventCut(){

}


