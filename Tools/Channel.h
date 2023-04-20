#ifndef CHANNEL_H
#define CHANNEL_H

#include <map>

#include "TString.h"
#include "EventCut.h"

class Channel {

public:

    Channel(TString name, TString latex,const TString cut);
    ~Channel();
 
    inline TString GetName(){return m_name;};
    inline TString GetLatex(){return m_latex;};
    inline EventCut *GetChannelCut(){return m_cut;};

    void AddCut(TString cut, TString name);
    inline std::vector<EventCut*> GetCut(){return m_vector_cut;};

    inline EventCut *GetExtraCut(){return m_extra;};
    void AddSampleCut(TString sample,TString cut,TString weight);
    inline TString GetSampleWeight(TString sample){return m_weight[sample];};   
    inline TString GetSampleCut(TString sample){return m_cut_sample[sample];};
    void AddVerticalLine(TString var,std::vector<float> values);
    void AddArrows(TString var,std::vector<float> values);
    void AddText(TString,std::vector< std::map<TString, float > > values);

    inline std::map<TString,std::vector<float> > GetVerticalLines(){return m_vertical_lines;};
    inline std::map<TString,std::vector<float> > GetArrows(){return m_arrows;};
    inline std::map<TString,std::vector< std::map<TString, float > > > GetText(){return m_text;};

    inline void SetTotalSystematic(float up,float down);
    inline std::vector<float> GetTotalSystematic(){return m_total_systematic;};


private:

    TString m_name;
    TString m_latex;
    std::vector<float> m_total_systematic;   
    std::vector<EventCut*> m_vector_cut;

    EventCut *m_cut;
    EventCut *m_extra;
    std::map<TString,TString> m_weight;
    std::map<TString,TString> m_cut_sample;    
    std::map<TString,std::vector<float> > m_vertical_lines;
    std::map<TString,std::vector<float> > m_arrows;
    std::map<TString,std::vector< std::map<TString, float > > > m_text;

};

#endif // CHANNEL_H
