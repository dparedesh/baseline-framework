#ifndef PHYSICSPROCESS_H
#define PHYSICSPROCESS_H


#include <vector>
#include <map>

#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"

#include "Yields.h"
#include "PhysicsSample.h"

class PhysicsProcess {

public:

    PhysicsProcess();
    ~PhysicsProcess();

    
    void AddSample(TString file,TString path,float scale_factor);

    void SetStyle(TString name,TString title,TString LabelLatex,int color,int line,TString process);
    inline void SetWeight(TString weight){m_weight=weight;};
    inline void SetCut(TString cut){m_cut=cut;};
    inline void SetTree(TString tree){m_tree=tree;};
    inline TString GetCut(){return m_cut;};
    inline TString GetWeight(){return m_weight;};
    inline TString GetTree(){return m_tree;};

    inline void AddCutFlow(std::vector<TH1F*> cutflow){m_cutflow=cutflow;};
    inline TH1F *GetCutFlow(int ch){return m_cutflow[ch];};
    
    inline void AddHistos(vector< vector<TH1F*> > histos){m_histos=histos;};
    inline TH1F* GetHistos(int fs,int var){return m_histos[fs][var];}
    inline std::vector< std::vector<TH1F*> > GetVectorHistos(){return m_histos;};    

    inline void AddHistos2D(vector< vector<TH2F*> > histos){m_histos_2D=histos;};
    inline TH2F* GetHistos2D(int fs,int var_2D){return m_histos_2D[fs][var_2D];}
    
    inline TString GetTitle(){return m_title;};
    inline TString GetName(){return m_name;};
    inline TString GetLatex(){return m_latex;};
    inline int GetColor(){return m_color;};
    inline int GetLine(){return m_line;};
    inline TString isDataSigBkg(){return m_process;};

    inline void SetYield(TString ch,Yields *yield){m_yield[ch]=yield;}
 
    inline double GetYield(TString ch){return m_yield[ch]->GetYield();}
    inline double GetStatistical(TString ch){return m_yield[ch]->GetStatistical();}

    inline std::vector<double> GetSystematic(TString ch){return m_yield[ch]->GetSystematic();}
    inline void AddPositionsMass(std::vector<TString> mass){m_positions_mass=mass;}
    inline std::vector<TString> GetPositionsMass(){return m_positions_mass;}
    inline void SetCombinedSignificance(double sig){m_comb_sig=sig;};
    inline double GetCombinedSignificance(){return m_comb_sig;};    
    inline void SetDictionary(int dict){m_dict=dict;};
    inline int GetDictionary(){return m_dict;};
    inline void SetIsDataDriven(bool isDD){m_is_datadriven=isDD;};
    inline bool GetIsDataDriven(){return m_is_datadriven;};
    void SetSystematicVariations(TString list_shapes,TString list_norm);
    inline std::map<TString,TString> GetSystematicsVariations(){return m_Systematics;};
    inline void AddSystematicVariations(std::map<TString,vector< vector<TH1F*> > > syst){m_histos_syst=syst;};
    inline TH1F* GetHistosSystematics(TString syst,int fs, int var){return m_histos_syst[syst][fs][var];};    
    inline std::map<TString,vector< vector<TH1F*> > > GetMapSystematics(){return m_histos_syst;};
    void SetSystematic(int fs,int var,TString ch,TString theory_unc);

    vector<PhysicsSample*> Samples;
 
private:

    std::vector<TString> FillVector(TString name);

    vector< vector<TH1F*> > m_histos;
    vector< vector<TH2F*> > m_histos_2D;
    std::vector<TH1F*> m_cutflow;
 
    std::map<TString,Yields*> m_yield;
    std::vector<TString> m_positions_mass;    
    std::map<TString,TString> m_Systematics;
    std::map<TString,TString> List_Systematics;
    std::map<TString,vector< vector<TH1F*> > > m_histos_syst;


    double m_comb_sig;
    TString m_title;
    TString m_name;
    TString m_latex;
    TString m_weight;
    TString m_cut;
    TString m_tree;
    int m_color;
    TString m_process;
    int m_line;
    int m_dict;
    bool m_is_datadriven;

};

#endif //PHYSICSPROCESS_H
