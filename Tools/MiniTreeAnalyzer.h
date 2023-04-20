#ifndef MINITREEANALYZER_H
#define MINITREEANALYZER_H

#include "string.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>

#include "Yields.h"
#include "EventCut.h"
#include "VariableDistr.h"
#include "VariableDistr2D.h"
#include "PhysicsProcess.h"
#include "Channel.h"

#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"


class MiniTreeAnalyzer {
 public:

    MiniTreeAnalyzer();
    ~MiniTreeAnalyzer();

    
    void AddChannel(Channel *p_ch);
    
    void AddCut(TString cut,TString name);
    //inline EventCut *GetFullCut(){return p_Cut;}
    void AddPreselectionCut(TString cut,TString name);
    inline EventCut *GetPreselectionCut(){return p_Cut_preselection;};
    
    inline void AddWeight(TString weight){m_weight=weight;}
    inline TString GetWeight(){return m_weight;}

    void SetTreeName(TString tree){m_tree=tree;};
    void SettingConstants();
    //inline void SetPrecision(int number){m_precision_table=number;};    

    void AddProcess(TString name,TString title,TString latex,int color,int line,TString process,bool isDD=false,TString cut="",TString weight="",TString tree="",int dict=0,TString syst_shape="",TString syst_norm="");
    void AddVariable(VariableDistr* var);
    void AddVariable2D(VariableDistr* var_x,VariableDistr* var_y);
    void PrintStatisticalUncertainty();


    void Execute();
    void MakeSkimming(TString cut,std::vector<TString> vars);
    void ComputeSystematicVariations();

    double ComputeSignificance(double sig,double bkg);
    double ComputeStatisticalSignificance(double sig,double bkg);
    double ComputeGeneralSignificance(double sig,double bkg,double ebkg);

    TGraphAsymmErrors* GetHistoPoisson(TH1F *pH);
    void GetATLAS(TString label,double x,double y,bool quad,double size);
    void GetLabel(double x,double y,TString label,float size);
    void SetATLASLabel(TString label="Internal",double x=0.2,double y=0.88,bool pos=false,double size=0.05); 
    void SetLuminosity(float lumi,TString label="",float unc=0);   
    TLegend* SetupLegend(float xl,float xu,float yl, float yu,int textfont,float textsize);
    inline void SetLegendParams(int ncolumns=1,TString total_title="Total SM",bool inSignal=true){m_ncolumns=ncolumns; m_total_leg=total_title; m_inSignal=inSignal;};

    void GetLabelMCTruthClassifier(TH1F* &pH,TString var);
    void GetLabelMCTruthClassifier(TH2F* &pH,TString var_x,TString var_y);
    std::vector<TString> ParticleMCTC(TString var);
    void DrawArrowsForPointsOutsideRange(TGraphAsymmErrors *pratio,double min,double max);
  
    inline vector<PhysicsProcess*> GetPhysicsProcesses(){return PhysProcess;};
    inline vector<Channel*> GetChannels(){return channels;};
    inline vector<VariableDistr*> GetVariables(){return variable;};
    //inline PhysicsProcess* GetPhysicsProcess(TString var,TString ch){}

    void SaveYieldsTables(bool save,TString folder="Tables");
    void SavePlots(bool save,TString folder="Plots");
    void SaveHistos(bool save,TString folder="Histos");
    void tokenizeString(std::string& str, char delim, std::vector<std::string>& tokens);
    void SetSystematicVariations(bool doSyst,TString list_shapes,TString list_norm,TString theory,TString OutputDir);
    std::map<TString,TString> GetSystematicsVariations(){return m_Systematics;};
    std::map<TString,TString> GetListSystematicsVariations(){return List_Systematics;};

    TH1F *CreateHisto(VariableDistr* var,TString ch);
    TH2F *CreateHisto2D(VariableDistr*  var_x,VariableDistr*  var_y, TString fs);
    inline void SetSignificance(double rel_syst=0.3,bool root_sig=true,bool stat_forum=false,bool standard=false,bool stat_standard=false){ m_rel_syst=rel_syst; m_root_sig=root_sig; m_stat_forum=stat_forum, m_standard=standard; m_stat_standard=stat_standard; doSignificance=true; };

    inline void SetDebugLevel(int debug=0){m_debug=debug;};

    inline void PlotNotStackedBkg(bool doShapes,std::map<TString,TString> norm,TString sample=""){m_doShapes=doShapes;m_norm_shapes=norm; m_sample_toNormalize_notStacked=sample;};

    TString datDir; // directory with data files
    TString sigDir; // directory with signal files
    TString bkgDir; // directory with background files
   
    bool doCutflow;

    bool doBkg;
    bool doSig;
    bool doData;
    bool doSkimming;
    bool printLog;  
    bool doRatioDataBkg;

    bool doSystematics;
    bool doSystematicVariations;
    bool doPlotDistributions;
    bool doStackSignal;
    bool doLogPlots;
    bool NormalizeToUnit;  //Each PhysicsProcess is scaled to 1 for plotting purposes
    


   vector<PhysicsProcess*> PhysProcess;
   vector<Channel*> channels;
   vector<VariableDistr*> variable;
   vector<VariableDistr2D*> variable2D;
  
   std::map<TString,Yields*> TotalBkg;
   std::map<TString,Yields*> TotalData;
   std::map<TString,std::map<TString,Yields*> > TotalSig;

   std::map< TString,std::map<TString,double> > Significance;
   std::map< TString,std::map<TString,double> > GeneralSignificance;
   std::map< TString,std::map<TString,double> > StatisticalSignificance;
   std::map< TString,std::map<TString,double> > RootSignificance;
   std::map< TString, double > CombinedSignificance;
   TH1F* GetSqrtTH1(TH1F *pH);
   TH1F* GetSquareTH1F(TH1F *pH);    


  private:

    
    std::vector<TString> FillVector(TString name);
    
    
    std::vector<TString> GetMassPoints(TString filename,int x, int y); 
    std::map<TString,TString> m_Systematics;
    std::map<TString,TString> List_Systematics;    
    std::map<TString,std::map<TString,TGraphAsymmErrors*> > m_total_bkg_syst;


    void LoopOverSamples();
    void SkimmingPerSample(TString file,TString cut,TString tree);   
    std::vector<TString> GetFullListSystematics(); 

    vector< vector<TH1F*> >   ReadInputFile(TString pFile,TString samplename,float sf,TString isdata,bool isDD,TString tree,TString process_weight,TString process_cut,int dict,std::vector<TH1F*>& cutflow);
    vector< vector<TH2F*> > ReadInputFile2D(TString pFile,TString samplename,float sf,TString isdata,bool isDD,TString tree,TString process_weight,TString process_cut,int dict);
    
    void ComputeYields();
    void GetSignificance();
    void PrintTableYields();
    void PrintInvertedTableYields();
    void PrintCutFlow();
    void BuildShapesNormalized();
    void BuildShapes();
    void MakePlotDistributions();
    void MakePlotDistributions2D();
    void SetHistoStyle(TH1F* &pH,int color,int line,TString process);
    void SetHistoStyle2D(TH2F* &pH,int color,int line,TString process);
    TString ConvertVariable(TString var,int dict);    
    TGraphAsymmErrors* GetTGraphTotalSystematic(std::vector<float> systematic, TH1F *pStat);
    double GetErrorRatio(float n, float d, float error_n, float error_d);
    TGraphAsymmErrors* GetRatio(TGraphAsymmErrors* num, TGraphAsymmErrors* den);
    TGraphAsymmErrors* GetTGraphErrorsSyst(TH1F *pStat,TH1F *up,TH1F *down);
    inline TGraphAsymmErrors* GetTotalBkgSystematics(TString ch,TString var){return m_total_bkg_syst[ch][var];};
    void SetTotalBkg();

    EventCut *p_Cut_preselection;
    vector<EventCut*> v_cuts;
    vector<TString> final_state;
    vector<TString> fs_latex;
    vector<TString> m_vars_skim;
    
    TString m_tree;
    TString m_weight;

    TString m_atlas_label;
    float m_atlas_x;
    float m_atlas_y;
    double m_atlas_size;
    bool m_atlas_pos;
    float m_luminosity;
    float m_luminosity_unc;
    TString m_luminosity_label;
    int m_precision_table;    
    int m_debug;

    bool m_save;
    TString m_histos_dir;
    bool m_save_plots;
    TString m_plots_dir;
    bool m_save_tables;
    TString m_tables_dir;
    TString OutputSystDir;
    TString file_theory_uncertainties;
    std::map<TString,TString> m_norm_shapes;


    int m_ncolumns;
    TString m_total_leg;
    bool m_inSignal;
    bool m_doShapes;
    bool m_sample_toNormalize_notStacked;

    double m_rel_syst;
    bool m_root_sig;
    bool m_stat_forum;
    bool m_standard;
    bool m_stat_standard;

    bool doSignificance;
};


#endif // MINITREEANALYZER_H
