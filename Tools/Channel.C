#include "Channel.h"

Channel::Channel(TString name, TString latex,const TString cut):
 m_name(name),
 m_latex(latex),
 m_weight(),
 m_total_systematic(),
 m_cut_sample(),
 m_vector_cut()

{

    m_cut=new EventCut(cut,name);
    m_extra=0;
    
    //std::cout << " Creating channel: '"<< name << "' with cut: " << cut << std::endl;
}

Channel::~Channel(){

    delete m_cut;
    delete m_extra;
  /*
   for (std::map<TString,TString>::iterator it=m_weight.begin(); it!=m_weight.end(); ++it){
        delete it->second;
    }

   for (std::map<TString,TString>::iterator it=m_cut_sample.begin(); it!=m_cut_sample.end(); ++it){
        delete it->second;
   }
   */

   for (unsigned int i=0; i<m_vector_cut.size(); i++) delete m_vector_cut[i];

   m_total_systematic.clear();
   m_weight.clear();
   m_cut_sample.clear();
   m_vertical_lines.clear();
   m_arrows.clear();
   m_text.clear();
   m_vector_cut.clear();

}
void Channel::AddCut(TString cut,TString name){
    
    //m_extra= new EventCut(cut,name);
    EventCut *local_cut= new EventCut(cut,name);

    m_vector_cut.push_back(local_cut);
 
    
    return;
}
void Channel::AddVerticalLine(TString var,std::vector<float> values){

 m_vertical_lines[var]=values; 


 return;
}
void Channel::AddArrows(TString var,std::vector<float> values){

  m_arrows[var]=values;

  return;
}
void Channel::AddText(TString var,std::vector< std::map<TString, float > > values){

  m_text[var]=values;

  return;
}
void Channel::SetTotalSystematic(float up, float down){

 m_total_systematic.push_back(up);
 m_total_systematic.push_back(down);

 return;
}
void Channel::AddSampleCut(TString sample,TString cut,TString weight){

   m_weight[sample]=weight;
   m_cut_sample[sample]=cut;


 return;
}
