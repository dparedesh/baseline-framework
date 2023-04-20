#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
//#include <iomanip>
//#include <vector>

#include "TPie.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TPaletteAxis.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TColor.h"

#include "Tools/Yields.C"
#include "Tools/Channel.C"
#include "Tools/EventCut.C"
#include "Tools/VariableDistr2D.C"
#include "Tools/VariableDistr.C"
#include "Tools/PhysicsSample.C"
#include "Tools/PhysicsProcess.C"

#include "Tools/MiniTreeAnalyzer.C"
#include "Tools/MyStyle.C"


void PrintOutput(Channel *Channels,vector<PhysicsProcess*> Processes,std::map<TString,Yields*> Bkg,std::map<TString,std::map<TString,double> > v_Significances);

void PlotsOffitial()
{    

     SetMyStyle();

     TString Rpc2L1b = "(nSigLep>=2 && nBJets20>=1 && nJets40>=6 && met/meff>0.25)";
     TString Rpc2L2b = "(nSigLep>=2 && nBJets20>=2 && nJets25>=6 && met>300000 && meff>1400000 && met/meff>0.14)";
     TString Rpc2L0b = "(nSigLep>=2 && nBJets20==0 && nJets40>=6 && met>200000 && meff>1000000 && met/meff>0.2)";
     TString Rpc3LSS1b = "(nSigLep>=3 && nBJets20>=1 && is3LSS>0 && is3LSSproc>0 && !isZee && met/meff>0.14)";
     TString Rpv2L0b = "(nSigLep>=2 && nBJets20>=0 && nJets40>=6 && meff>2600000)";
   
     TString VRWZ4j = "(nSigLep==3 && NlepBL==3 && nBJets20==0 && nJets25>=4 && mSFOS>81000 && mSFOS<101000 && met>50000 && met<250000 && meff>600000 && meff<1500000 ) && !"+Rpc2L1b+" && !"+Rpc2L2b+" && !"+Rpc2L0b+" && !"+Rpv2L0b;
     TString VRWZ5j = "(nSigLep==3 && NlepBL==3 && nBJets20==0 && nJets25>=5 && meff>400000 && met>50000 && mSFOS>81000 && mSFOS<101000 && meff<1500000 && met<250000) && !"+Rpc2L1b+" && !"+Rpc2L2b+" && !"+Rpc2L0b+" && !"+Rpv2L0b;
     TString VRttV = "(nSigLep>=2 && NlepBL>=2 && nBJets20>=1 && nJets40>=3 && meff>600000 && meff<1500000 && met<250000 && isSS30 && dRl1j>1.1 && SumBJetPt/SumJetPt>0.4 && met/meff>0.1) && !"+Rpc2L1b+" && !"+Rpc2L2b+" && !"+Rpc2L0b+" && !"+Rpv2L0b;


      Channel *ch_Rpc2L0b= new Channel("Rpc2L0b","Rpc2L0b",Rpc2L0b);
      Channel *ch_Rpc2L1b= new Channel("Rpc2L1b","Rpc2L1b",Rpc2L1b);
      Channel *ch_Rpc2L2b= new Channel("Rpc2L2b","Rpc2L2b",Rpc2L2b);
      Channel *ch_Rpc3LSS1b= new Channel("Rpc3LSS1b","Rpc3LSS1b",Rpc3LSS1b);
      Channel *ch_Rpv2L0b= new Channel("Rpv2L","Rpv2L",Rpv2L0b);

      Channel *ch_VRttV= new Channel("VRttV","VRttV",VRttV);
      Channel *ch_VRWZ4j= new Channel("VRWZ4j","VRWZ4j",VRWZ4j);
      Channel *ch_VRWZ5j= new Channel("VRWZ5j","VRWZ5j",VRWZ5j);

      MiniTreeAnalyzer testanalyzer;

       
      VariableDistr *var = new VariableDistr("nSigLep","nSigLep","nSigLep","Events",10,-10,10,true); 
      VariableDistr *var_njets = new VariableDistr("nJets25","nJets25","n_{jets}","Events",9,3.5,12.5,true);
      VariableDistr *v_meff = new VariableDistr("meff","meff","m_{eff} [GeV]","Events/",11,400,1500,true);

      MiniTreeAnalyzer analyzer;
      analyzer.datDir="/afs/cern.ch/user/d/dparedes/WorkCERN/Analysis_SUSY_SS3L_2018/HistFitterComplexFinal/SS3Lep/SS3L_HF/FullRun2/prepare/"; 
      analyzer.bkgDir="/eos/user/d/dparedes/SUSYComplex/";
      analyzer.sigDir="/eos/user/d/dparedes/SUSYComplex/";
      analyzer.SetTreeName("evtel");
      analyzer.SetATLASLabel("Internal"); 
      analyzer.SetLuminosity(139,"139 fb^{-1}");   
      analyzer.SetDebugLevel(1); 

       
      //analyzer.AddChannel(ch_Rpc2L0b);
      analyzer.AddChannel(ch_Rpc2L1b);
      //analyzer.AddChannel(ch_Rpc2L2b);
      //analyzer.AddChannel(ch_Rpc3LSS1b);
      //analyzer.AddChannel(ch_Rpv2L0b);
      
      //analyzer.AddChannel(ch_VRttV);
      //analyzer.AddChannel(ch_VRWZ4j);
     

      //analyzer.AddChannel(ch_VRWZ5j);

      //analyzer.AddVariable(var);
      analyzer.AddVariable(var_njets);
      //analyzer.AddVariable(v_meff);

      
      //analyzer.AddProcess("InputList/Bkg.txt","Multiboson","Multiboson",28,1,"isBkg",false,"","","Multiboson_nom",0);
      //analyzer.AddProcess("InputList/Bkg.txt","ttV","t#bar{t}V",kOrange,1,"isBkg",false,"","","ttV_nom",0);
      analyzer.AddProcess("InputList/Bkg.txt","OtherMultiboson","OtherMultiboson",kGreen+2,1,"isBkg",false,"","","OtherMultiboson_nom",3);
      analyzer.AddProcess("InputList/Bkg.txt","ttH","ttH",kOrange,1,"isBkg",false,"","","ttH_nom",3); 
      analyzer.AddProcess("InputList/Bkg.txt","RareTop_NottH","RareTop (no ttH)",kCyan+1,1,"isBkg",false,"","","RareTop_NottH_nom",3);
      analyzer.AddProcess("InputList/Bkg.txt","ttW","t#bar{t}W",kViolet-9,1,"isBkg",false,"","","ttW_nom",3);
      analyzer.AddProcess("InputList/Bkg.txt","ttZ","t#bar{t}Z",kAzure+7,1,"isBkg",false,"","","ttZ_nom",3);
      analyzer.AddProcess("InputList/Bkg.txt","WZ","WZ",kOrange+1,1,"isBkg",false,"","","WZ_nom",3);
      analyzer.AddProcess("InputList/Other.txt","FakesMC","FakesMC",19,1,"isBkg",false,"","","FakesMC_nom",3);         

      analyzer.AddProcess("InputList/Data.txt","Data","Data",1,1,"isData",false,"","","data_nom",3);

      analyzer.AddProcess("InputList/signal.txt","Signal","Signal",kRed,2,"isSig",false,"","","",3);


     
      analyzer.AddWeight("totweight*lumiScaling"); //("wmu_nom*wel_nom*wtrig_nom*wjet_nom*mcweight*wpu_nom_bkg*wpu_nom_sig*MC_campaign_weight*lumiScaling");  
      //analyzer.SetPrecision(2);
      analyzer.doPlotDistributions=true;
      analyzer.SetSignificance(0.3,true,true);
      analyzer.SavePlots(true,"Plots");
      analyzer.doLogPlots=false;
      analyzer.doRatioDataBkg=true;
      analyzer.SaveHistos(true,"Histos");
      analyzer.SaveYieldsTables(true,"Tables");
      //analyzer.SetSystematicVariations(false,"Syst_shape_listWZ.txt","Syst_norm_short.txt","TheoryUnc.root","OutputSystematics"); //"Syst_norm_list.txt","OutputSystShapes");

      analyzer.doCutflow=false;
      analyzer.Execute();
 
    

    //=====================================================================================
    //Getting out values:
    vector<PhysicsProcess*> Processes=analyzer.GetPhysicsProcesses();
    vector<Channel*> Channels=analyzer.GetChannels();
    std::map<TString,Yields*> Bkg=analyzer.TotalBkg;
    std::map< TString,std::map<TString,double> > v_Significances=analyzer.RootSignificance;


    for (unsigned int j=0; j<Channels.size(); j++){
        std::cout << "--- In channel : " << Channels[j]->GetName() << std::endl;
        PrintOutput(Channels[j],Processes,Bkg,v_Significances);
    }
 



return;
}
void PrintOutput(Channel *Channels,vector<PhysicsProcess*> Processes,std::map<TString,Yields*> Bkg,std::map<TString,std::map<TString,double> > v_Significances){


  std::cout << "---- Total bkg: " << Bkg[Channels->GetName()]->GetYield() << " +/- " << Bkg[Channels->GetName()]->GetStatistical()   << std::endl;
  for (unsigned int i=0; i<Processes.size(); i++){

     std::cout << "-- Process :" << Processes[i]->GetTitle() << ",     yield: "<<  Processes[i]->GetYield(Channels->GetName()) << "+/-" << Processes[i]->GetStatistical(Channels->GetName())   << std::endl;

     if (Processes[i]->isDataSigBkg()=="isSig"){
        std::cout << "..... Significance: " << v_Significances[Processes[i]->GetName()][Channels->GetName()] << std::endl;
     }
  }



  return;
}


