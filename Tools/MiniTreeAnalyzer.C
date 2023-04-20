#include "MiniTreeAnalyzer.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <regex>

#include "TFile.h"
#include "TTree.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TArrow.h"
#include "TText.h"
#include <TROOT.h>
#include <TStyle.h>
#include "RooStats/NumberCountingUtils.h"

#include "RooHistError.h"

//#include "boost/algorithm/string.hpp"

using namespace RooStats; 

MiniTreeAnalyzer::MiniTreeAnalyzer():
datDir(""),
sigDir(""),
bkgDir(""),

OutputSystDir(""),
file_theory_uncertainties(""),

doBkg(false),
doData(false),
doSig(false),
doCutflow(false),
doPlotDistributions(false),
doStackSignal(false),
doSkimming(false),
doLogPlots(false),
doRatioDataBkg(false),
doSystematics(false),
doSystematicVariations(false),
NormalizeToUnit(false),
printLog(false),


m_save(false),
m_histos_dir(""),
m_save_tables(false),
m_tables_dir(""),
m_save_plots(false),
m_plots_dir(""),
m_precision_table(3),

doSignificance(false),
m_tree(""),
m_weight(""),

Significance(),
GeneralSignificance(),
StatisticalSignificance(),
RootSignificance(),
CombinedSignificance(),
m_Systematics(),
List_Systematics(),

m_atlas_label(""),
m_atlas_x(0),
m_atlas_y(0),
m_atlas_pos(false),
m_atlas_size(0),
m_luminosity(0),
m_luminosity_unc(0),
m_luminosity_label(""),

m_ncolumns(1),
m_total_leg("Total SM"),
m_inSignal(true),

final_state(),
fs_latex(),

m_norm_shapes(),

m_sample_toNormalize_notStacked(""),
m_doShapes(false),

TotalBkg(),
TotalData(),
TotalSig(),
m_total_bkg_syst(),
m_vars_skim(),

variable(),
variable2D(),
channels(),
v_cuts(),
//p_Cut(),
p_Cut_preselection(),
PhysProcess()

{}

MiniTreeAnalyzer::~MiniTreeAnalyzer()
{
    
    if (p_Cut_preselection) delete p_Cut_preselection;
    
    for (unsigned int i=0; i<PhysProcess.size(); i++){
        delete PhysProcess[i];
    }
 
    /*for (unsigned int i=0; i<channels.size(); i++){
        delete channels[i];
    }*/
    for (unsigned int i=0; i<v_cuts.size(); i++){
        delete v_cuts[i];
    }
    for (unsigned int i=0; i<variable.size(); i++){
        delete variable[i];
    }
    for (unsigned int i=0; i<variable2D.size(); i++){
        delete variable2D[i];
    }
    
    
    for (std::map<TString,Yields*>::iterator it=TotalBkg.begin(); it!=TotalBkg.end(); ++it){
        delete it->second;
    }


    for (std::map<TString, std::map<TString,TGraphAsymmErrors*> >::iterator it=m_total_bkg_syst.begin(); it!=m_total_bkg_syst.end(); ++it){  
 
            for (std::map<TString,TGraphAsymmErrors*>::iterator it_s=it->second.begin(); it_s != it->second.end(); ++it_s){
         
                delete it_s->second;

             }
    }


   /* for (std::map<TString,Yields*>::iterator it=TotalData.begin(); it!=TotalData.end(); ++it){
        delete it->second;
    }*/
    
    
    /*for (std::map<TString, std::map<TString,Yields*> >::iterator it=TotalSig.begin(); it!=TotalSig.end(); ++it){
        
	 delete it->second;
	
	 for (std::map<TString,Yields*>::iterator it_second=it->second.begin(); it_second != it->second.end(); ++it_second){
         
             delete it_second->second;

         }
    }*/
    
   
    
    PhysProcess.clear();
    channels.clear();
    v_cuts.clear();
    variable.clear();
    variable2D.clear();
    final_state.clear();
    fs_latex.clear();
    
    
    
    TotalBkg.clear();
    TotalData.clear();
    TotalSig.clear();
    m_total_bkg_syst.clear();   
 
    Significance.clear();
    GeneralSignificance.clear();
    StatisticalSignificance.clear();
    RootSignificance.clear(); 
    CombinedSignificance.clear();
    m_Systematics.clear();
    List_Systematics.clear();

}


void MiniTreeAnalyzer::MakeSkimming(TString cut,std::vector<TString> vars){


 m_vars_skim=vars;
 
 std::cout<< "######## Starting skimming with cut: "<< cut << std::endl;

 for (unsigned int i=0; i<PhysProcess.size(); i++){

        if (m_debug>-1) cout << "### STARTING PROCESS:"<< PhysProcess[i]->GetName() << endl;
        
        for (unsigned int j=0; j<PhysProcess[i]->Samples.size(); j++){

        	SkimmingPerSample(PhysProcess[i]->Samples[j]->GetFile(),cut,PhysProcess[i]->GetTree());


                 if (doSystematicVariations){

		   std::cout << "--- In systematics #########" << std::endl;	

          	   std::map<TString,TString> systematics_list=PhysProcess[i]->GetSystematicsVariations();

		   if (systematics_list.empty() && !PhysProcess[i]->GetIsDataDriven() && PhysProcess[i]->isDataSigBkg()!="isData") systematics_list=GetSystematicsVariations();

                   for (std::map<TString,TString>::iterator it=systematics_list.begin(); it!=systematics_list.end(); ++it){

                        TString local_tree=PhysProcess[i]->GetTree();
                        TString local_weight=PhysProcess[i]->GetWeight();

 			std::cout << "------- Systematics name: " << local_tree << std::endl;

                        if (it->second=="isShape") local_tree.ReplaceAll("nom",it->first); //PhysProcess[i]->GetTitle()+"_"+it->first;
                        else if (it->second=="isNorm") continue;

			std::cout << "------- Systematics name: " << local_tree << std::endl;

			SkimmingPerSample(PhysProcess[i]->Samples[j]->GetFile(),cut,local_tree);


                   } //end for syst 

            } // end if doSyst




        } // end of samples

 }//end for PhysProcess    



 return;
}

void MiniTreeAnalyzer::SkimmingPerSample(TString pFile,TString cut,TString tree){

    TFile *pInput = TFile::Open(pFile);
    if (!pInput){
        printf("-- File %s is missing\n", pInput->GetName());
        exit(1);
    }
    else {
        if (m_debug>0) cout << "-- Opening file :" << pFile << endl;
    }


    std::vector<std::string> temp;
  
   
    std::string sample(pFile.Data());
   
    tokenizeString(sample,'/',temp); 

    TString ntree=tree;
    if (tree.IsNull()) ntree=m_tree;


    TTree* local_Tree = (TTree*)pInput->Get(ntree);
    local_Tree->SetBranchStatus("*",kFALSE);

    for (unsigned int i=0; i<m_vars_skim.size(); i++) local_Tree->SetBranchStatus(m_vars_skim[i],kTRUE);


    if (m_debug>0) std::cout<< "-- Tree entries:" << local_Tree->GetEntries() << std::endl;

    TH1F *pHfile;

    if (!pFile.Contains("data") && !pFile.Contains("Fakes") && !pFile.Contains("ChargeFlip")) pHfile = (TH1F*)pInput->Get("hCutFlow"); 

    TString newfile=bkgDir+"/Skimming/"+temp[temp.size()-1];
    std::cout <<"--- Creatting file : " << newfile << std::endl; 

    TFile *pOutput = new TFile(newfile,"RECREATE");

  
    TTree* new_tree = local_Tree->CopyTree(cut.Data());


    if (pHfile && !pFile.Contains("data") && !pFile.Contains("Fakes")) pHfile->Write();

    std::cout <<"## entries in skimmed tree : "<< new_tree->GetEntries() << std::endl; 

    new_tree->Print();
    new_tree->Write();
    pOutput->Close();
    pInput->Close();

    delete pOutput;
    delete pInput;


 return;
}


void MiniTreeAnalyzer::Execute()
{

 //std::cout << "************* Starting execution to read MiniTrees--" << std::endl;

    SettingConstants();
 
    LoopOverSamples();

    //SetTotalBkg();


    if (doSystematicVariations) ComputeSystematicVariations();    

    ComputeYields();
  
    //PrintStatisticalUncertainty();

    if (doCutflow) PrintCutFlow();

    
    if (doPlotDistributions) {
         if (NormalizeToUnit) BuildShapesNormalized();
         else if (m_doShapes) BuildShapes();
         else MakePlotDistributions();

         if (variable2D.size()>0) MakePlotDistributions2D();
    }
    

  
  if (m_save_tables){    
    for (unsigned int i=0; i<variable.size(); i++){
    
        if (variable[i]->GetYield()){
            PrintTableYields();
            PrintInvertedTableYields();
            break;
        }
    }
  }
 

 //std::cout<<"-- End Execute() " << std::endl;

 if (m_debug>0) std::cout << "### Ending execution of MiniTreeAnalyser..." << std::endl;
 if (m_debug>0) std::cout <<"#######################################################################" << std::endl;

 return;
}


void MiniTreeAnalyzer::BuildShapes(){


 std::cout <<"######## DEBUG --> WARNING: When providing the samples in AddSample(),  the denominator samples must be provided first than numerator. Otherwise you will probably get a segmentation fault!  #####" << std::endl;

 TString PlotsDir=m_plots_dir;

 if (m_save_plots) system("mkdir -p "+PlotsDir);

 

 for (unsigned int fs=0; fs<final_state.size(); fs++){ //fs=final_state.size()-1
        for (unsigned int var=0; var< variable.size(); var++){


            TCanvas *pC = new TCanvas(PlotsDir+"_h_"+variable[var]->GetTitle()+"_"+final_state[fs],"h_"+variable[var]->GetTitle()+"_"+final_state[fs],700,600);

            TPad* c_plot;
            TPad* c_ratio;

            c_plot = new TPad("plot","",0.,0.3,1.0,1.0);
            c_plot->SetRightMargin(0.08);
            c_plot->SetLeftMargin(0.15);
            c_plot->SetBottomMargin(0.0);
            c_plot->SetTopMargin(0.05);
            c_plot->Draw();

            c_ratio = new TPad("ratio","",0.,0,1,0.3);
            c_ratio->SetRightMargin(0.08);
            c_ratio->SetLeftMargin(0.15);
            c_ratio->SetBottomMargin(0.35);

            c_ratio->SetTopMargin(0.0);
            c_ratio->Draw();
            c_ratio->SetGridy(kTRUE);

            c_plot->cd();

            float h=4.5*0.0775;
            
            TLegend *pLeg=SetupLegend(0.5,0.9-h,0.5+0.35,0.9,42,0.045); //0.57
            pLeg->SetNColumns(m_ncolumns);

            TString kind_of_process="";

            if (variable[var]->plotLogX()) c_plot->SetLogx();

            float max=0;

            //TH1F *pNormalize=0
            //TH1F *hratio=0;
            std::map<TString,TH1F*> hist_den;
            std::map<TString,TH1F*> hist_num;
            std::vector<TH1F*> hist_ratios;

            for (unsigned process=0; process< PhysProcess.size(); process++){

                TH1F *h_den=0;
                TH1F *h_num=0;

                TH1F *pH=(TH1F*)PhysProcess[process]->GetHistos(fs,var)->Clone();
                kind_of_process=PhysProcess[process]->isDataSigBkg();

                pH->SetLineColor(PhysProcess[process]->GetColor());
                pH->SetLineWidth(2);
                pH->SetLineStyle(PhysProcess[process]->GetLine());

                //int bin_max=pH->GetMaximum
                //pH->SetMaximum(1);

                pH->GetXaxis()->SetLabelSize(0.00);
                pH->GetXaxis()->SetTitleSize(0.00);
                pH->GetXaxis()->SetTitleOffset(1.0);

                pH->GetYaxis()->SetTitleOffset(1.1);
                pH->GetYaxis()->SetTitleSize(0.07);
                pH->GetYaxis()->SetLabelSize(0.06);

                pH->SetMinimum(0);


                if (process==0){
                  int bin_max=pH->GetMaximumBin();
                  pH->SetMaximum(pH->GetBinContent(bin_max)*1.5);
                  pH->Draw("E");
                  pH->Draw("hist,same");

                }
                else {
                    pH->Draw("hist,same");
                    pH->Draw("E,same");

                }

                pLeg->AddEntry(pH,PhysProcess[process]->GetLatex(),"l");

                TString process_name=PhysProcess[process]->GetTitle();

                std::cout <<" ######DEBUG:" << process_name << std::endl;


                if (m_norm_shapes.count(process_name)>0){
                        h_den=(TH1F*)pH->Clone();
                        hist_den[process_name]=h_den;
                        std::cout <<"### debug: den :" << process_name << ": " << h_den->Integral() << std::endl;

                }
                else {
                       h_num=(TH1F*)pH->Clone();
                       hist_num[process_name]=h_num;
                       std::cout <<"### debug: num :" << process_name << ": " << h_num->Integral() << std::endl;

                }

                /*if (process==0) pNormalize=(TH1F*)pH->Clone();
                else {
                    std::cout << "### DEBUG : COMPUTING RATIO "<< std::endl;
                    hratio=(TH1F*)pH->Clone();
                    hratio->Divide(pNormalize);
                }*/

            }  // end physics process

            pLeg->Draw("same");


            //GetATLAS(m_atlas_label,m_atlas_x,m_atlas_y,m_atlas_pos,m_atlas_size);
            //GetLabel(0.2,0.7,fs_latex[fs],0.055); //0.2,0.78,0.04
            GetATLAS(m_atlas_label,0.2,0.86,m_atlas_pos,0.06);
            GetLabel(0.20,0.79,"#sqrt{s} = 13 TeV, "+m_luminosity_label,0.06);
            GetLabel(0.2,0.72,fs_latex[fs],0.053);


            gPad->RedrawAxis();

            pC->Update();


            c_ratio->cd();

            if (variable[var]->plotLogX()) c_ratio->SetLogx();

            
            TH1F *pLine = (TH1F*)hist_den.begin()->second->Clone();
            for (int i=1; i<=pLine->GetNbinsX(); i++) pLine->SetBinContent(i,1);
            pLine->SetLineWidth(2);
            pLine->SetLineColor(1);

            pLine->SetMarkerSize(1.2); //0.8
            pLine->SetLineWidth(2);
            //pLine->GetYaxis()->SetRangeUser(0.05,0.25);

            pLine->GetXaxis()->SetTitleOffset(1.0);
            pLine->GetYaxis()->SetTitleOffset(0.5);
            pLine->GetYaxis()->SetNdivisions(505); //404
            pLine->GetXaxis()->SetTitleSize(0.15);
            pLine->GetYaxis()->SetTitleSize(0.14);
            pLine->GetXaxis()->SetLabelSize(0.14);
            pLine->GetYaxis()->SetLabelSize(0.14);
            pLine->GetXaxis()->SetTickSize(0.07);
            pLine->GetYaxis()->SetTitle("Ratio");
            pLine->GetYaxis()->SetRangeUser(0.,0.5);
            c_ratio->SetGridy();
            pLine->Draw("histo");
            
            for ( auto it = m_norm_shapes.begin(); it != m_norm_shapes.end(); ++it  ){
               std::cout <<" @@ DEBUG : "<< it->first << ", " << it->second << std::endl;

               TH1F *pH_local=(TH1F*)hist_num[it->second]->Clone();
               pH_local->Divide(hist_den[it->first]);

               pH_local->SetFillColor(1);
               pH_local->SetFillStyle(3004);
               pH_local->Draw("pl,same");

            }

            /*
            hratio->SetFillColor(1);
            hratio->SetFillStyle(3004);
            hratio->Draw("pl,same");*/
            gPad->RedrawAxis();

            if (m_save_plots) pC->SaveAs(PlotsDir+"/all_"+variable[var]->GetTitle()+"_"+final_state[fs]+".eps");


         } // end of var
 }// end of final_State   


return;
}



void MiniTreeAnalyzer::BuildShapesNormalized(){

 TString PlotsDir=m_plots_dir;

 if (m_save_plots) system("mkdir -p "+PlotsDir);


 for (unsigned int fs=0; fs<final_state.size(); fs++){ //fs=final_state.size()-1
        for (unsigned int var=0; var< variable.size(); var++){
           

            TCanvas *pC = new TCanvas(PlotsDir+"_h_"+variable[var]->GetTitle()+"_"+final_state[fs],"h_"+variable[var]->GetTitle()+"_"+final_state[fs],700,600);
            pC->cd();

            
            TLegend *pLeg=new TLegend(0.61,0.68,0.63+0.28,0.93);  //0.55;
            pLeg->SetNColumns(m_ncolumns);
            
           
            TString kind_of_process="";
            
            
            for (unsigned process=0; process< PhysProcess.size(); process++){
                
                TH1F *pH=(TH1F*)PhysProcess[process]->GetHistos(fs,var)->Clone();
                kind_of_process=PhysProcess[process]->isDataSigBkg();
                
                //SetHistoStyle(pH,PhysProcess[process]->GetColor(),PhysProcess[process]->GetLine(),kind_of_process);
   
   		pH->SetLineColor(PhysProcess[process]->GetColor());
		pH->SetLineWidth(2);
                pH->SetLineStyle(PhysProcess[process]->GetLine());  

                if (pH->Integral()!=0) pH->Scale(1.0/pH->Integral());

               
                int bin_max=pH->GetMaximumBin();
                pH->SetMaximum(pH->GetBinContent(bin_max)*1.3);     


      		if (process==0) pH->Draw("hist");
		else pH->Draw("hist,same");
		
		pLeg->AddEntry(pH,PhysProcess[process]->GetLatex(),"l");
		
            }  // end physics process
    
    		pLeg->Draw("same");
    
 
           GetATLAS(m_atlas_label,0.2,0.88,m_atlas_pos,0.041);
           GetLabel(0.2,0.83,"#sqrt{s} = 13 TeV, "+m_luminosity_label,0.041);
           GetLabel(0.2,0.78,fs_latex[fs],0.038);

    
            if (m_save_plots) pC->SaveAs(PlotsDir+"/all_"+variable[var]->GetTitle()+"_"+final_state[fs]+".eps");
    

         } // end of var
}// end of final_State	 



 return;
}
void MiniTreeAnalyzer::MakePlotDistributions(){

 if (m_debug>0) std::cout <<"### MakePlotDistributions()" << std::endl;

 bool transpose=false;

 TString PlotsDir=m_plots_dir;

 if (m_save_plots) system("mkdir -p "+PlotsDir);


 TString HistosDir=m_histos_dir;

 if (m_save) system("mkdir -p "+HistosDir);
 


    for (unsigned int fs=0; fs<final_state.size(); fs++){ //fs=final_state.size()-1

        std::map<TString,std::vector<float> >v_lines=channels[fs]->GetVerticalLines();
        std::map<TString,std::vector<float> >v_arrows=channels[fs]->GetArrows();
        std::map<TString,std::vector< std::map<TString, float > > > v_text=channels[fs]->GetText();


        for (unsigned int var=0; var< variable.size(); var++){
           
            THStack *p_bkg= new THStack("bkg_"+final_state[fs]+"_"+variable[var]->GetTitle(),"");
            TH1F *data_histo=0;
            
            TCanvas *pC = new TCanvas(PlotsDir+"_h_"+variable[var]->GetTitle()+"_"+final_state[fs],"h_"+variable[var]->GetTitle()+"_"+final_state[fs],700,600);
          
	    TPad* c_plot;
            TPad* c_ratio;


           if (doRatioDataBkg && doData){

      		c_plot = new TPad("plot","",0.,0.3,1.0,1.0);
                c_plot->SetRightMargin(0.08);
                c_plot->SetLeftMargin(0.15);
                c_plot->SetBottomMargin(0.0);
                c_plot->SetTopMargin(0.05);
 		c_plot->Draw();

      		c_ratio = new TPad("ratio","",0.,0,1,0.3);
                c_ratio->SetRightMargin(0.08);
                c_ratio->SetLeftMargin(0.15);
                c_ratio->SetBottomMargin(0.35);

      		c_ratio->SetTopMargin(0.0); 
      		c_ratio->Draw();
      		c_ratio->SetGridy(kTRUE);
	   }	
           else {
                c_plot = new TPad("plot","",0.,0.,1.0,1.0);
                c_plot->SetRightMargin(0.08);
                c_plot->SetLeftMargin(0.15);
                c_plot->SetBottomMargin(0.15);
                c_plot->SetTopMargin(0.05);
                c_plot->Draw();
 
		c_ratio = 0;
	   }

           c_plot->cd();



            TLegend *pLeg;

            //setupLegend(nleg, 2, 0.53, 0.90, 0.35, -1, 0.07, fontsize=0.045)
            float h=4.5*0.0775;

            if (doSig) pLeg=new TLegend(0.51,0.58,0.63+0.28,0.93);
            if (doRatioDataBkg && doData) pLeg=SetupLegend(0.53,0.9-h,0.53+0.35,0.9,42,0.045); //0.57
            else pLeg=new TLegend(0.61,0.68,0.63+0.28,0.93);  //0.55;

            TLegend *pLegTrans;


            if (doRatioDataBkg) pLegTrans = SetupLegend(0.53,0.9-h,0.53+0.35,0.9,42,0.045);
            else pLegTrans=SetupLegend(0.50,0.9-h,0.50+0.40,0.9,42,0.038);

            pLeg->SetNColumns(m_ncolumns);
            pLegTrans->SetNColumns(m_ncolumns);
           
            TString kind_of_process="";
            
            float maxSignal=0;
            float maxTemp=0;
            
           
            std::vector<TH1F*> vec_signal;
            
            int bkg_counter=0;

            for (unsigned process=0; process< PhysProcess.size(); process++){
                
                TH1F *pH=(TH1F*)PhysProcess[process]->GetHistos(fs,var)->Clone();
                kind_of_process=PhysProcess[process]->isDataSigBkg();
                
                SetHistoStyle(pH,PhysProcess[process]->GetColor(),PhysProcess[process]->GetLine(),kind_of_process);
   
               
                
                if (kind_of_process=="isBkg" && doBkg==true) {
                      p_bkg->Add((TH1F*)pH->Clone());
		      if (PhysProcess[process]->GetLatex() != "TRAMPA") bkg_counter=bkg_counter+1;
                }
                else if (kind_of_process=="isData") data_histo=(TH1F*)pH->Clone();
                else if (kind_of_process=="isSig"){
                    maxTemp=pH->GetMaximum();
                    if (maxTemp>maxSignal) maxSignal=maxTemp;
                    
                    vec_signal.push_back((TH1F*)pH->Clone());
                }

            }  // end physics process
    
            
            TH1F* Temporal=0;
            TH1F* h_base=0;
            TH1F* Tempo1=0;
            
            if (doBkg) Temporal=(TH1F*)p_bkg->GetStack()->Last();
            else Temporal=(TH1F*)PhysProcess[0]->GetHistos(fs,var)->Clone();
            
            
            Tempo1= (TH1F*)Temporal->Clone();
            h_base=Tempo1;
            
            h_base->SetLineColor(0);
            h_base->SetFillColor(0);
            h_base->SetMarkerColor(0);
            
            
            float maxStack=0;
            
            if (doStackSignal){
                std::vector<THStack*> signal_stack;
                
                for (unsigned int i=0; i<vec_signal.size(); i++){
                    
                    THStack *phstack= new THStack( "Temp_"+final_state[fs]+"_"+variable[var]->GetTitle(),"");
                    phstack->Add((TH1F*)h_base->Clone());
                    phstack->Add((TH1F*)vec_signal[i]->Clone());
                    
                    if (maxStack < phstack->GetMaximum()) maxStack=phstack->GetMaximum();
                    
                    
                }
            }
            
            
            
            float max=h_base->GetMaximum()+h_base->GetBinError(h_base->GetMaximumBin());
            
            if (doData && (data_histo->GetMaximum()+data_histo->GetBinError(data_histo->GetMaximumBin())> max)) max=data_histo->GetMaximum()+data_histo->GetBinError(data_histo->GetMaximumBin());
            
            
            if (m_debug>0) std::cout << "-- Plotting variable :" << variable[var]->GetName() << std::endl;
            
            if (maxSignal>max) max=maxSignal;
            
            h_base->SetMaximum(2*max);
            h_base->SetMinimum(0.0001);
            
            if (variable[var]->GetName().Contains("truthType") || variable[var]->GetName().Contains("truthOrig")) GetLabelMCTruthClassifier(h_base,variable[var]->GetName()); 

            /*
             h_base->GetXaxis()->SetTitleOffset(1.2);
             h_base->GetYaxis()->SetTitleOffset(1.4);
             h_base->GetXaxis()->SetTitleSize(0.00);
             h_base->GetYaxis()->SetTitleSize(0.05);
             h_base->GetXaxis()->SetLabelSize(0.00);
             h_base->GetYaxis()->SetLabelSize(0.05);
             */

            if (h_base->GetXaxis()->GetBinLowEdge(h_base->GetNbinsX()+1)>1500) h_base->GetXaxis()->SetNdivisions(505);

            h_base->Draw("hist"); //DrawCopy
          
            if (doStackSignal){
                std::vector<THStack*> signal_stack;
                
		int counter=0;

                for (unsigned int i=0; i<vec_signal.size(); i++){
                    
		  if (counter<4){
                    THStack *phstack= new THStack( "_"+final_state[fs]+"_"+variable[var]->GetTitle(),"");
                    phstack->Add((TH1F*)h_base->Clone());
                    phstack->Add((TH1F*)vec_signal[i]->Clone());
                    
                    phstack->Draw("hist,same");
                
		  }	

		  counter=counter+1;
    
                }
            }
            
            
            
            
            TGraphAsymmErrors *pSystematic=0;
            
            if (doBkg){
                p_bkg->Draw("hist,same");
                h_base->SetLineWidth(2);
                h_base->SetLineColor(1);
		h_base->DrawCopy("hist,same");                


                h_base->SetFillStyle(3004); //3254
                h_base->SetFillColor(1);
                h_base->DrawCopy("same E2");

                
                if (doSystematicVariations) {
              
                   //if (doSystematics) pSystematic=GetTGraphTotalSystematic(channels[fs]->GetTotalSystematic(),h_base);
                    pSystematic=GetTotalBkgSystematics(channels[fs]->GetName(),variable[var]->GetName()); 
    
                    pSystematic->SetFillColor(1);
                    pSystematic->SetFillStyle(3004);
                   // pSystematic->Draw("e2,same");
                }
                
            }
        
	    if (pSystematic !=0 ){
                 TGraphAsymmErrors *pSystematicCopy = (TGraphAsymmErrors*)pSystematic->Clone();  
                 pSystematicCopy->SetFillColor(1);
                 pSystematicCopy->SetFillStyle(3004);
                 pSystematicCopy->Draw("same e2");

             }


            
            double binning=h_base->GetBinWidth(1);
            //TString str_bin=ConvertInt(binning);
            //TString name_file=HistosDir+"/"+"Histo_"+variable[var]->GetName()+"_"+final_state[fs]+"_"+CutLevel+"_bin"+str_bin+".root";
            //if (variable[var]->IsArray())name_file=HistosDir+"/"+"Histo_"+variable[var]->GetName()+"_"+final_state[fs]+"_"+CutLevel+".root";
            
            //system("rm -rf "+name_file);
            
            
            
            //TFile *pFile = new TFile(name_file,"RECREATE");
            
            for (unsigned int process=0; process<PhysProcess.size(); process++){
                
                kind_of_process=PhysProcess[process]->isDataSigBkg();
                TH1F *pH=(TH1F*)PhysProcess[process]->GetHistos(fs,var)->Clone();
                SetHistoStyle(pH,PhysProcess[process]->GetColor(),PhysProcess[process]->GetLine(),kind_of_process);
                
                if (kind_of_process=="isData") {
                    pLeg->AddEntry(pH,PhysProcess[process]->GetLatex(),"p");
                    pLegTrans->AddEntry(pH,PhysProcess[process]->GetLatex(),"p");
                }
            }

	     if (doSystematicVariations) {
                     pSystematic->SetLineWidth(2);
                     pLeg->AddEntry(pSystematic,m_total_leg,"lf");  
          
             }      
             else {
                     h_base->SetLineWidth(2);
                     pLeg->AddEntry(h_base,m_total_leg,"lf");
             }



            TString name_file=HistosDir+"/Histo1D_"+variable[var]->GetTitle()+"_"+final_state[fs]+".root";
            TFile *pFile;

            if (m_save) pFile = new TFile(name_file,"RECREATE");


            std::vector<PhysicsProcess*> bkg_trans(bkg_counter);

            

            int local_pos=0;

            for (signed int process=PhysProcess.size()-1; process>-1; process--){
                
                kind_of_process=PhysProcess[process]->isDataSigBkg();
                TH1F *pH=(TH1F*)PhysProcess[process]->GetHistos(fs,var)->Clone();
                pH->SetDirectory(0);
                SetHistoStyle(pH,PhysProcess[process]->GetColor(),PhysProcess[process]->GetLine(),kind_of_process);
                
                
                if (kind_of_process=="isBkg" && PhysProcess[process]->GetLatex() != "TRAMPA") {
                           pLeg->AddEntry(pH,PhysProcess[process]->GetLatex(),"f");

                        

                          //int local_pos=process;
                       /*
                           if (PhysProcess[process]->GetTitle()=="ZZ") local_pos=0;
                           else if (PhysProcess[process]->GetTitle()=="ttV") local_pos=1;
                           else if (PhysProcess[process]->GetTitle()=="Fakes") local_pos=2;
                           else if (PhysProcess[process]->GetTitle()=="ChargeFlip") local_pos=3;
                           else if (PhysProcess[process]->GetTitle()=="WZ") local_pos=4;
                           else if (PhysProcess[process]->GetTitle()=="Rare") local_pos=5;
                           else if (PhysProcess[process]->GetTitle()=="WW") local_pos=6;
                        */

			  if (m_debug>1) std::cout << "############## DEBUG:: LOCAL POSITION:" << PhysProcess[process]->GetTitle() << " --> " << local_pos << std::endl;
                         

                          bkg_trans[local_pos]=PhysProcess[process];

                          local_pos+=1;

                 }


                 if (m_save){
                    pH->SetName(PhysProcess[process]->GetTitle());
                    pH->Write();
                  }

       
            } // end of PhysProcess
            


            for (unsigned int k=0; k<bkg_trans.size(); k++){


                TH1F *pH=(TH1F*)bkg_trans[k]->GetHistos(fs,var)->Clone();
                pH->SetDirectory(0);
                SetHistoStyle(pH,bkg_trans[k]->GetColor(),bkg_trans[k]->GetLine(),bkg_trans[k]->isDataSigBkg());


                 if (k==0){

                       pLegTrans->AddEntry(pH,bkg_trans[k]->GetLatex(),"f");

                       if (doSystematicVariations) pLegTrans->AddEntry(pSystematic,m_total_leg,"lf");
                       else pLegTrans->AddEntry(h_base,m_total_leg,"lf");
                 }
                 else pLegTrans->AddEntry(pH,bkg_trans[k]->GetLatex(),"f");



            }
            




            if (m_save) pFile->Close();            
            

	    int counter_sig=0;



            TLegend *plegsignal;
            if (doRatioDataBkg) plegsignal=SetupLegend(0.525,0.46,0.52+0.35,0.56,42,0.045);
            else  plegsignal = SetupLegend(0.49,0.46,0.49+0.4,0.56,42,0.037); //SetupLegend(0.50,0.9-h,0.50+0.50,0.9,42,0.038);/

            for (unsigned int process=0; process<PhysProcess.size(); process++){
                
                kind_of_process=PhysProcess[process]->isDataSigBkg();
                TH1F *pH=(TH1F*)PhysProcess[process]->GetHistos(fs,var)->Clone();
                SetHistoStyle(pH,PhysProcess[process]->GetColor(),PhysProcess[process]->GetLine(),kind_of_process);
                
                if (kind_of_process=="isSig") { 
             
		  if (counter_sig<4){

                    if (m_inSignal){
                       pLeg->AddEntry(pH,PhysProcess[process]->GetLatex(),"lf");
                       pLegTrans->AddEntry(pH,PhysProcess[process]->GetLatex(),"lf");
                    }

		    plegsignal->AddEntry(pH,PhysProcess[process]->GetLatex(),"lf");
       
                    if(!doStackSignal) pH->Draw("hist,same");
                    
		  }

		  counter_sig=counter_sig+1;

                } // end if isSig   
            }
            
            
 	    TGraphAsymmErrors *pData=0;
            
            if (doData) {


                pData=GetHistoPoisson(data_histo);
                pData->SetMarkerStyle(20);
                pData->Draw("p,same");
                pData->SetLineWidth(1);
                pData->SetMarkerSize(1.1);
                pData->SetName("Data");         


		if (doRatioDataBkg){
      		   /* h_base->GetYaxis()->SetLabelSize(0.06);   // 0.1
                    h_base->GetYaxis()->SetTitleSize(0.06);   // 0.12
                    h_base->GetYaxis()->SetTitleOffset(0.8);*/

                   h_base->GetXaxis()->SetLabelSize(0.00);
		   h_base->GetXaxis()->SetTitleSize(0.00);
   	           h_base->GetXaxis()->SetTitleOffset(1.0);

                   h_base->GetYaxis()->SetTitleOffset(1.1);
                   h_base->GetYaxis()->SetTitleSize(0.07);
                   h_base->GetYaxis()->SetLabelSize(0.06);


		}
                else{

                   h_base->GetXaxis()->SetLabelSize(0.04);
                   h_base->GetXaxis()->SetTitleSize(0.05);
                   h_base->GetXaxis()->SetTitleOffset(1.1);

                   h_base->GetYaxis()->SetTitleOffset(1.1);
                   h_base->GetYaxis()->SetTitleSize(0.05);
                   h_base->GetYaxis()->SetLabelSize(0.04);
                   h_base->GetYaxis()->SetNdivisions(510);

		}
                
        
            }// end doData
    
            
 
            float factor=1.3;

            if (v_lines.count(variable[var]->GetTitle())==1){

		std::vector<float> lines=v_lines[variable[var]->GetTitle()];

		for (unsigned int li=0; li<lines.size(); li++){
	       
			TLine *pline=new TLine(lines[li],0.001,lines[li],factor*max);
			pline->SetLineStyle(1);
       			pline->SetLineWidth(3);
                        pline->SetLineColor(12);

			pline->Draw("same");
   
                        TLine *pclone=(TLine*)pline->Clone();
                        pclone->Draw("same");
		}
	    }
            if (v_arrows.count(variable[var]->GetTitle())==1){

                std::vector<float> arrows=v_arrows[variable[var]->GetTitle()];           

                 for (unsigned int ar=0; ar<arrows.size()-1; ar+=2){

			TArrow *parrow= new TArrow(arrows[ar],factor*max,arrows[ar+1],factor*max,0.01); //,0.05,"|>");
                        parrow->SetLineStyle(1);
                        parrow->SetLineWidth(3);
                        parrow->SetLineColor(12);

                        parrow->Draw();

                 }
            }
           
            if (v_text.count(variable[var]->GetTitle())==1){

                std::vector< std::map<TString, float > > text=v_text[variable[var]->GetTitle()];

                for (unsigned int tx=0; tx<text.size(); tx++){

			for (std::map<TString,float>::iterator it=text[tx].begin(); it!=text[tx].end(); ++it){

			    TLatex *pText = new TLatex(); 
                            pText->SetTextColor(12);

                            if (doRatioDataBkg && doData) pText->SetTextSize(0.04);
			    else pText->SetTextSize(0.026);
                            pText->DrawLatex(it->second,(factor+0.05)*max,it->first.Data());
                        }

                 }

            } 


            pLeg->SetTextFont(42);
            pLegTrans->SetTextFont(42); 
           
            if (!transpose) pLeg->Draw("same");
            else pLegTrans->Draw("same");

            if (!m_inSignal){
              plegsignal->SetTextFont(42);
              plegsignal->Draw("same");     
            }

            if (doRatioDataBkg && doData){
                 GetATLAS(m_atlas_label,0.2,0.86,m_atlas_pos,0.06); 
	         GetLabel(0.20,0.79,"#sqrt{s} = 13 TeV, "+m_luminosity_label,0.06);
                 GetLabel(0.2,0.72,fs_latex[fs],0.053);  
            }
            else {           
           	 GetATLAS(m_atlas_label,0.2,0.88,m_atlas_pos,0.041);
          
            	 GetLabel(0.2,0.83,"#sqrt{s} = 13 TeV, "+m_luminosity_label,0.041);
            	 GetLabel(0.2,0.78,fs_latex[fs],0.038);
            }    
            
	    gPad->RedrawAxis();
	    

            pC->Update();

            if (doRatioDataBkg && doData){

		TH1F *hratio = (TH1F*)data_histo->Clone();

		TH1F *pLine = (TH1F*)h_base->Clone();


		for (int i=1; i<=pLine->GetNbinsX(); i++) pLine->SetBinContent(i,1);
		pLine->SetLineWidth(2);

	
               

                c_ratio->cd();

 		hratio->Divide(h_base);
                hratio->SetMarkerSize(1.1); //0.8
                hratio->SetLineWidth(1);
                hratio->GetYaxis()->SetRangeUser(0.0001,1.999);
                //hratio->GetYaxis()->CenterTitle();

	/*	hratio->GetXaxis()->SetLabelSize(0.14);  // 0.12
       		hratio->GetXaxis()->SetTitleSize(0.17); // 0.15
       		hratio->GetXaxis()->SetTitleOffset(0.9); //0.08 
       		hratio->GetYaxis()->SetLabelSize(0.13);   // 0.1
       		hratio->GetYaxis()->SetTitleSize(0.13);   // 0.12
       		hratio->GetYaxis()->SetTitleOffset(0.33);*/

                hratio->GetXaxis()->SetTitleOffset(1.0);
                hratio->GetYaxis()->SetTitleOffset(0.5);
                hratio->GetYaxis()->SetNdivisions(404);
                hratio->GetXaxis()->SetTitleSize(0.15);
                hratio->GetYaxis()->SetTitleSize(0.14);
                hratio->GetXaxis()->SetLabelSize(0.14);
                hratio->GetYaxis()->SetLabelSize(0.14);
                hratio->GetXaxis()->SetTickSize(0.07);




                hratio->GetYaxis()->SetTitle("Data/SM");	
 		hratio->SetLineColor(0);
	        hratio->SetMarkerColor(0);
                //hratio->GetYaxis()->SetNdivisions(405);
                hratio->Draw("p");
		
                  
		  
		TGraphAsymmErrors *pErrorBand=0;
                TGraphAsymmErrors *pratio=0;

                if (doSystematicVariations){
                           pErrorBand=GetRatio(pSystematic,pSystematic);
                           pratio=GetRatio(pData,pSystematic);
                } 
                else {

                   TGraphAsymmErrors *pStatTemp=new TGraphAsymmErrors(h_base);

		   pErrorBand=GetRatio(pStatTemp,pStatTemp);
		   pratio=GetRatio(pData,pStatTemp);
                }
   
		pErrorBand->SetFillColor(1);
                pErrorBand->SetFillStyle(3004);
                pErrorBand->Draw("same,e2");	
                   
	        pratio->SetMarkerStyle(20);
                pratio->SetMarkerSize(1.1);
                pratio->SetLineWidth(1);
                pratio->Draw("p,same");


                Double_t* xcoor=pErrorBand->GetX();
                int npoints=pErrorBand->GetN();

                DrawArrowsForPointsOutsideRange(pratio,xcoor[0],xcoor[npoints-1]);


		TLine* line1 = new TLine(hratio->GetXaxis()->GetXmin(),1.,hratio->GetXaxis()->GetXmax(),1.);
		line1->SetLineWidth(1);
  		line1->Draw("same");             


          }


           

            if (m_save_plots){
                  pC->SaveAs(PlotsDir+"/all_"+variable[var]->GetTitle()+"_"+final_state[fs]+".eps");
                  pC->SaveAs(PlotsDir+"/all_"+variable[var]->GetTitle()+"_"+final_state[fs]+".png");
                  pC->SaveAs(PlotsDir+"/all_"+variable[var]->GetTitle()+"_"+final_state[fs]+".pdf");
                  TFile *pCanvas = new TFile(PlotsDir+"/"+variable[var]->GetTitle()+"_"+final_state[fs]+".root","RECREATE");
                  pC->Write();
                  pCanvas->Close(); 
           }

	    
            pC->Update();  	
    
	    /// Make log plots here:
	
          if (doLogPlots){
     
	    float max_log=10000*h_base->GetMaximum()+h_base->GetBinError(h_base->GetMaximumBin());
	    
	    h_base->SetMaximum(max_log);
	    h_base->SetMinimum(0.001);


	    if (v_lines.count(variable[var]->GetTitle())==1){

                std::vector<float> lines=v_lines[variable[var]->GetTitle()];

                for (unsigned int li=0; li<lines.size(); li++){

                        TLine *pline=new TLine(lines[li],0.001,lines[li],0.001*max_log);
                        pline->SetLineStyle(2);
                        pline->SetLineWidth(2);

                        pline->Draw("same");
                }
            }	


	   	    
	    c_plot->SetLogy();
	    pC->Update();

	    if (m_save_plots) pC->SaveAs(PlotsDir+"/all_"+variable[var]->GetTitle()+"_"+final_state[fs]+"_log.eps");
  
                
         }    	    
	    
	    
        } // end var
    }  // end fs
    

    if (m_debug>0) std::cout << "#### Ending MakePlotDistributions() "<< std::endl;    
    
    return;
}

TLegend* MiniTreeAnalyzer::SetupLegend(float xl,float yl,float xu, float yu,int textfont,float textsize){

 TLegend *pLeg= new TLegend(xl,yl,xu,yu);

 pLeg->SetTextFont(textfont);
 pLeg->SetTextSize(textsize);

 return pLeg;
}
void MiniTreeAnalyzer::DrawArrowsForPointsOutsideRange(TGraphAsymmErrors* plot,double xmin,double xmax){

 Double_t* y=plot->GetY();
 Double_t* x=plot->GetX();
 int n=plot->GetN();  


  int counter=0;
  for (int j=0; j<n; j++){

       if ( y[j] >= 2 || y[j] <0 ){
             if (m_debug>1) std::cout << "################## DEBUGGING .....  " << j+1 	<< std::endl;
             float yl=1.74, yu=1.995;

             if (y[j]<0){
                 yl=0.26;
		 yu=0.05;
             }
             TArrow *parrow = new TArrow(x[j],yl,x[j],yu,0.01);
             parrow->SetLineWidth(2);
             //parrow->SetLineColor(12);
             if (counter==0) parrow->Draw();  
             else parrow->Draw("same");

             float range=xmax-xmin;

             float x_pos=x[j]+0.005*range;
             float y_text=int(y[j]*10.0)/10.0;

             string text=to_string(y_text);
             TString newtext="";

             int temp=0;
             for (int i=0; i<text.size(); i++){
  
                char chr=text[i];

                 if (temp<=1){
                    if (chr=='.' || temp>0) {
                            temp=temp+1;
                            newtext.Append(text[i]);
		    }
                    else newtext.Append(text[i]);	
                }
                else break;

             }
             if (m_debug>1) std::cout << "-- val : " << newtext << std::endl;

             TText *ptext = new TText(x_pos,yl,newtext);
             ptext->SetTextFont(62);
             ptext->SetTextSize(0.07);
             ptext->Draw();
     } 

  }



  return;
}
TGraphAsymmErrors* MiniTreeAnalyzer::GetRatio(TGraphAsymmErrors* num, TGraphAsymmErrors* den){

 // This function return the ratio of num and den where the uncertainties on the denominator are set to 0.

 TGraphAsymmErrors *ratio = new TGraphAsymmErrors;


 int n_den=den->GetN();
 int n_num=num->GetN();

 Double_t* y_num=num->GetY();
 Double_t* x_num=num->GetX();


Double_t* x_den=den->GetX();
Double_t* y_den=den->GetY();


 int counter=0;

 for (int j=0; j<n_den; j++){

    for (int k=0; k<n_num; k++){	

         if (x_den[j]==x_num[k]){
  
	    float x=x_num[k];
     	    float x_error=den->GetErrorXhigh(k);

     	    float y=y_num[k]/y_den[j];

     	    float num_up=num->GetErrorYhigh(k);
     	    float den_up=0; //den->GetErrorYhigh(j);

     	    float num_dn=num->GetErrorYlow(k);
     	    float den_dn=0; //den->GetErrorYlow(j);

     	    float y_up=GetErrorRatio(y_num[k],y_den[j],num_up,den_up);
     	    float y_down=GetErrorRatio(y_num[k],y_den[j],num_dn,den_dn);	

            if (strncmp(num->GetName(),"Data",4)==0) x_error=0;

	    ratio->SetPoint(counter,x,y);
	    ratio->SetPointError(counter,x_error,x_error,y_down,y_up);

	    //std::cout <<" num y: "<< y_num[k] << ", num y up: " << num_up << ", num_ y down:" << num_dn << std::endl;

	    //std::cout << " x , y, dx, dyup, dydn: "<< x << ", " << y << ", " << x_error << "," << y_up << ", " << y_down << std::endl;

	     counter=counter+1;
	     break;
	 }

    }	

 }



 return ratio;
}



double MiniTreeAnalyzer::GetErrorRatio(float n, float d, float error_n, float error_d){

 double error2=0;

 error2=pow(error_n/d,2)+pow(n/(d*d),2)*pow(error_d,2);


 return TMath::Sqrt(error2);
}
TGraphAsymmErrors* MiniTreeAnalyzer::GetTGraphTotalSystematic(std::vector<float> systematic, TH1F *pStat){

 TGraphAsymmErrors *pTotalSystematicBkg = new TGraphAsymmErrors;

 float up=systematic[0];
 float down=systematic[1];

 TH1F *pHUpTotal = (TH1F*)pStat->Clone();
 TH1F *pHDownTotal = (TH1F*)pStat->Clone();
  
 pHUpTotal->Reset();
 pHDownTotal->Reset();



 for (signed int bin=0; bin < pStat->GetNbinsX(); bin++){

      double stat2=pow(pStat->GetBinError(bin+1),2);

      double y=pStat->GetBinContent(bin+1);
      double x=pStat->GetBinCenter(bin+1);
      double dx=0.5*pStat->GetBinWidth(bin+1);

      double syst2_up=pow(up*y,2);
      double syst2_down=pow(down*y,2);

      double total_up=TMath::Sqrt(stat2+syst2_up);
      double total_down=TMath::Sqrt(stat2+syst2_down);

      pHUpTotal->SetBinContent(bin+1,TMath::Sqrt(syst2_up));
      pHDownTotal->SetBinContent(bin+1,TMath::Sqrt(syst2_down));

      pTotalSystematicBkg->SetPoint(bin,x,y);
      pTotalSystematicBkg->SetPointError(bin,dx,dx,total_down,total_up);

 }

 std::cout << "---- Total Bkg Systematic -> "<< pHDownTotal->Integral() << " / +" << pHUpTotal->Integral() << std::endl; 

  
 

  return pTotalSystematicBkg;

}
void MiniTreeAnalyzer::MakePlotDistributions2D(){

  TString HistosDir=m_histos_dir;

  if (m_save){
   system("mkdir -p "+HistosDir);
  }

 TString PlotsDir=m_plots_dir;

 if (m_save_plots){
   system("mkdir -p "+PlotsDir);
 }



  for (unsigned int fs=0; fs<final_state.size(); fs++){ //fs=final_state.size()-1
   for (unsigned int var_2D=0; var_2D< variable2D.size(); var_2D++){ 


       if (m_debug>0) std::cout << "### In Final state :" << final_state[fs] << std::endl; 

       TString name_file=HistosDir+"/"+"Histo2D_"+variable2D[var_2D]->GetVarY()->GetTitle()+"_vs_"+variable2D[var_2D]->GetVarX()->GetTitle()+"_"+final_state[fs]+"_.root";
       TFile *pFile;
       
       if (m_save) pFile = new TFile(name_file,"RECREATE");


       //TCanvas *pC = new TCanvas(PlotsDir+"h2D_"+variable2D[var_2D]->GetVarY()->GetTitle()+"_vs_"+variable2D[var_2D]->GetVarX()->GetTitle()+"_"+final_state[fs],"h2D_"+variable2D[var_2D]->GetVarY()->GetTitle()+"_vs_"+variable2D[var_2D]->GetVarX()->GetTitle()+"_"+final_state[fs],5);
       //pC->cd();

       //TLegend *pLeg= new TLegend(0.63,0.70,0.63+0.28,0.93); 

       TString var_x=variable2D[var_2D]->GetVarX()->GetName();
       TString var_y=variable2D[var_2D]->GetVarY()->GetName();

       for (unsigned process=0; process< PhysProcess.size(); process++){
           //TLegend *pLeg= new TLegend(0.63,0.70,0.63+0.28,0.93);

           TCanvas *pC = new TCanvas(PlotsDir+"h2D_"+variable2D[var_2D]->GetVarY()->GetTitle()+"_vs_"+variable2D[var_2D]->GetVarX()->GetTitle()+"_"+final_state[fs]+"_"+PhysProcess[process]->GetTitle(),"h2D_"+variable2D[var_2D]->GetVarY()->GetTitle()+"_vs_"+variable2D[var_2D]->GetVarX()->GetTitle()+"_"+final_state[fs]+"_"+PhysProcess[process]->GetTitle(),5);    
           pC->cd();
           pC->SetRightMargin(0.15);


	   TH2F *pH=(TH2F*)PhysProcess[process]->GetHistos2D(fs,var_2D)->Clone();
	   pH->SetDirectory(0);
       
           if (m_debug>0) std::cout << "-- MakePlotDistributions2D() : " << PhysProcess[process]->GetName() << "  -->  " << pH->Integral() << std::endl;


           if (var_x.Contains("truthType") || var_y.Contains("truthType") || var_x.Contains("truthOrig") || var_y.Contains("truthOrig") ) GetLabelMCTruthClassifier(pH,var_x,var_y);

           TString kind_of_process=PhysProcess[process]->isDataSigBkg();

	   SetHistoStyle2D(pH,PhysProcess[process]->GetColor(),PhysProcess[process]->GetLine(),kind_of_process);	   

            gStyle->SetPaintTextFormat("4.2f");
            gStyle->SetTextSize(0.001);
            gStyle->SetTextFont(42);
            pH->SetMarkerColor(kRed+1);
            pH->SetMarkerSize(0.8);

            pH->GetZaxis()->SetTitleOffset(1);

            if (pH->Integral()==1)pH->GetZaxis()->SetTitle("[%]");

            pH->GetZaxis()->SetNdivisions(505);


           if (pH->GetYaxis()->GetBinLowEdge(pH->GetNbinsY()+1)>1500) pH->GetYaxis()->SetNdivisions(505);
           if (pH->GetXaxis()->GetBinLowEdge(pH->GetNbinsX()+1)>1500) pH->GetXaxis()->SetNdivisions(505);


	   /*if (process==0) pH->DrawCopy("BOX,TEXT45"); //BOX
           else {
	      if (kind_of_process=="isdata") pH->DrawCopy("SCAT, same");	
	      else pH->DrawCopy("BOX,TEXT45"); //BOX,same
	   }*/
 
           /*if (process==0) pH->DrawCopy("SCAT"); //BOX
           else pH->DrawCopy("SCAT,SAME"); //BOX,same 
           */ 
           pH->DrawCopy("COLZ,TEXT45");

	   
	   //pLeg->AddEntry(pH,PhysProcess[process]->GetLatex(),"lfp"); //lf

          if (m_save){
	    pH->SetName(PhysProcess[process]->GetTitle());
            pH->Write();
          }

          GetATLAS(m_atlas_label,0.2,0.88,m_atlas_pos,0.040);
          GetLabel(0.2,0.83,"#sqrt{s} = 13 TeV, "+m_luminosity_label,0.040);
          GetLabel(0.2,0.78,PhysProcess[process]->GetLatex(),0.038);


          //pLeg->Draw("same");
          gPad->RedrawAxis();

          if (m_save_plots) pC->SaveAs(PlotsDir+"/all_"+variable2D[var_2D]->GetVarY()->GetTitle()+"_vs_"+variable2D[var_2D]->GetVarX()->GetTitle()+"_"+final_state[fs]+"_"+PhysProcess[process]->GetTitle()+".eps");


       }//end for PhysProcess

       


     // if (m_save_plots) pC->SaveAs(PlotsDir+"/all_"+variable2D[var_2D]->GetVarY()->GetTitle()+"_vs_"+variable2D[var_2D]->GetVarX()->GetTitle()+"_"+final_state[fs]+".eps");


      if (m_save) pFile->Close();
   } //end for var_2D

 }//end for fs


 return;
}




void MiniTreeAnalyzer::SetLuminosity(float lumi,TString label="",float unc=0){

    m_luminosity=lumi;
    m_luminosity_unc=unc;
    m_luminosity_label=label;

    if (m_luminosity_label=="") m_luminosity_label=std::to_string(lumi)+" fb^{-1}";

 return;
}
void  MiniTreeAnalyzer::SetATLASLabel(TString label="Internal",double x=0.2,double y=0.88,bool pos=false,double size=0.05){

 m_atlas_label=label;
 m_atlas_x=x;
 m_atlas_y=y;
 m_atlas_pos=pos;
 m_atlas_size=size;

 return;
}
void MiniTreeAnalyzer::GetATLAS(TString label,double x,double y,bool quad,double size){
    
    TLatex *pLat1 = new TLatex;
    
    pLat1->SetNDC();
    pLat1->SetTextColor(1);
    pLat1->SetTextSize(size); //0.05
    pLat1->SetTextFont(72);
 
  /*   
    float sep=0.19; //0.38
    float alt=0.85; //0.87
    float alts=0.07;
    */
    
    pLat1->DrawLatex(x,y,"ATLAS");
    pLat1->SetTextFont(42);
    //pLat1->SetTextSize(0.04);

    if (size==0.05) pLat1->SetTextSize(0.04);
    
    TString full_label="#font[42]{"+label+"}";
    
    double xx=0.17;
    if (!quad) xx=0.12;
    
    pLat1->DrawLatex(x+xx,y,full_label);
        
    
    
    return;
}
void MiniTreeAnalyzer::GetLabel(double x,double y,TString label,float size){
    
    TLatex *pLat1 = new TLatex;
    
    pLat1->SetNDC();
    pLat1->SetTextColor(1);
    
    pLat1->SetTextFont(42);
    pLat1->SetTextSize(size);
    
    pLat1->DrawLatex(x,y,label);
    
    
    return;
}

TGraphAsymmErrors* MiniTreeAnalyzer::GetHistoPoisson(TH1F *pH){
    
    TGraphAsymmErrors *pAsym = new TGraphAsymmErrors;
    
    int counter=0;

    for (signed int i=0; i<pH->GetNbinsX(); i++ ){
        
        float x=pH->GetBinCenter(i+1);
        
        Int_t y=pH->GetBinContent(i+1);
        
        if (y>0){
            
            Double_t yl=0, yu=0;
            
            RooHistError::instance().getPoissonInterval(y,yl,yu,1);
            
            pAsym->SetPoint(counter,x,y);
            pAsym->SetPointError(counter,0,0,y-yl,yu-y);
            counter=counter+1;
        }
        
    }
    
    
    return pAsym;
}


void MiniTreeAnalyzer::SetHistoStyle(TH1F* &pH,int color,int line,TString process)
{
    
    if (process=="isBkg"){
        pH->SetFillColor(color);
        pH->SetMarkerStyle(0);
        pH->SetLineColor(1);
        //pH->SetLineColor(color+1);
        //pH->SetLineWidth(2);
    }
    else if (process=="isData"){
        pH->SetMarkerColor(color);
        pH->SetLineColor(color);
        pH->SetMarkerStyle(20);
        pH->SetMarkerSize(1.2);
        pH->SetLineWidth(2);
    }
    
    else if (process=="isSig"){
        pH->SetMarkerColor(color);
        pH->SetLineColor(color);
        if (doStackSignal){
            pH->SetFillStyle(3002);
            pH->SetFillColor(color);
        }
        pH->SetLineWidth(5);
        pH->SetLineStyle(line);
    }
    
    
    return;
}

void MiniTreeAnalyzer::SetHistoStyle2D(TH2F* &pH,int color,int line,TString process)
{

  if (process != "isData"){
    pH->SetFillColor(color);
    pH->SetMarkerColor(color);
    //pH->SetMarkerStyle(0);
    //pH->SetLineColor(1);
    pH->SetLineWidth(0);
  }
  else{
   // pH->SetMarkerColor(color);
    pH->SetFillColor(0);
    pH->SetLineColor(color);
    pH->SetMarkerColor(color);
   // pH->SetMarkerStyle(20);
   // pH->SetMarkerSize(1.1);
    //pH->SetLineWidth(1);
  }



return;
}


void MiniTreeAnalyzer::PrintTableYields(){
      
  TString output=m_tables_dir;

  system("rm -rf "+output);
  system("mkdir -p "+output);
    
  std::ofstream out;
  TString filename=output+"/Table_Yields.tex";
  out.open(filename,std::ofstream::out);
    
    out << "\\begin{table}[htdp]"  << std::endl;
    out << "\\begin{center}"  << std::endl;
    out << "\\begin{tabular}{l";
    cout.flush();
	for (unsigned int fs=0; fs<channels.size(); fs++){
            if (fs==channels.size()-1) out << "|c}" << std::endl;
            else out << "|c";
        } 

    out << "\\hline" << std::endl;
    out << "\\hline" << std::endl;
    cout.flush();


    out << "\\multirow{2}{*}{Process}  & \\multicolumn{"<< channels.size()<<"}{c}{Channel} \\\\ " << std::endl;
    out << "\\cline{2-" << channels.size()+1 << "}" << std::endl;

    for (unsigned int fs=0; fs<channels.size(); fs++){

	TString local=channels[fs]->GetLatex().ReplaceAll("#","\\");
        if (local.Contains("\\")) local="$"+local+"$";

        if (fs==channels.size()-1) out << " & " << local << " \\\\"   << std::endl;
        else out << " & "<< local;
    }	 

    out << "\\hline" << std::endl;

    for (unsigned int process=0; process< PhysProcess.size(); process++){
        
        if (PhysProcess[process]->isDataSigBkg() != "isSig" && PhysProcess[process]->isDataSigBkg() !="isData" ){
            
 	    TString local=PhysProcess[process]->GetLatex().ReplaceAll("#","\\");
            if (local.Contains("\\")) local="$"+local+"$";            

	     out << local;
		
              for (unsigned int fs=0; fs<final_state.size(); fs++){
                
                   float yield=PhysProcess[process]->GetYield(channels[fs]->GetName());
                   float unc=PhysProcess[process]->GetStatistical(channels[fs]->GetName());		               
                   std::vector<double> syst = PhysProcess[process]->GetSystematic(channels[fs]->GetName());             


		   if (fs==final_state.size()-1) {
                        if (doSystematicVariations)  out << " & $"<< std::fixed << std::setprecision(m_precision_table) << yield << "\\pm "  
								  << std::fixed << std::setprecision(m_precision_table) << unc << "^{+"<< std::fixed << std::setprecision(m_precision_table) << syst[0] << "}_{-"  
								  << std::fixed << std::setprecision(m_precision_table) << syst[1] << "}" << "$ \\\\" << std::endl;
                        else  out << " & $"<< std::fixed << std::setprecision(m_precision_table) << yield << "\\pm "  << std::fixed << std::setprecision(m_precision_table) << unc  << "$ \\\\" << std::endl;
                   }

                   else {

                        if (doSystematicVariations) out << " & $"<< std::fixed << std::setprecision(m_precision_table) << yield << "\\pm "  
								 << std::fixed << std::setprecision(m_precision_table) << unc << "^{+"<< std::fixed << std::setprecision(m_precision_table) << syst[0] << "}_{-"
								 << std::fixed << std::setprecision(m_precision_table) << syst[1] << "}" << "$";
                        else out << " & $"<< std::fixed << std::setprecision(m_precision_table) << yield << "\\pm "  << std::fixed << std::setprecision(m_precision_table) << unc << "$";	       	        
                
		   }

              } // end fs

        } //end if 
    }// end process bkg and data
    
   out << "\\hline" << std::endl;
   
   out << " Total ";
   for (unsigned int fs=0; fs<final_state.size(); fs++){

                   float yield=TotalBkg[channels[fs]->GetName()]->GetYield();
                   float unc=TotalBkg[channels[fs]->GetName()]->GetStatistical();
 		   std::vector<double> syst=TotalBkg[channels[fs]->GetName()]->GetSystematic();

		   if (fs==final_state.size()-1) {
                     if (doSystematicVariations)  out << " & $"<< std::fixed << std::setprecision(m_precision_table) << yield << "\\pm "
                                                                  << std::fixed << std::setprecision(m_precision_table) << unc << "^{+"<< std::fixed << std::setprecision(m_precision_table) << syst[0] << "}_{-"
                                                                  << std::fixed << std::setprecision(m_precision_table) << syst[1] << "}" << "$ \\\\" << std::endl;
                     else out << " & $"<< std::fixed << std::setprecision(m_precision_table) << yield << "\\pm "  << std::fixed << std::setprecision(m_precision_table) << unc  << "$ \\\\" << std::endl;
                   }

                   else {

		    if (doSystematicVariations) out << " & $"<< std::fixed << std::setprecision(m_precision_table) << yield << "\\pm "  
                                                                 << std::fixed << std::setprecision(m_precision_table) << unc << "^{+"<< std::fixed << std::setprecision(m_precision_table) << syst[0] << "}_{-"
                                                                 << std::fixed << std::setprecision(m_precision_table) << syst[1] << "}" << "$";	
                     else out << " & $"<< std::fixed << std::setprecision(m_precision_table) << yield << "\\pm "  << std::fixed << std::setprecision(m_precision_table) << unc << "$";	
                   }
   }

   if (doData){
   	out << "\\hline" << std::endl;
	out << " Data ";
	for (unsigned int fs=0; fs<final_state.size(); fs++){

                   float yield=TotalData[channels[fs]->GetName()]->GetYield();
                   float unc=TotalData[channels[fs]->GetName()]->GetStatistical();


                   if (fs==final_state.size()-1) {
                         out << " & $" << int(yield) << "$ \\\\" << std::endl;
                   }

                   else out << " & $" << int(yield) << "$";

        }
   }

    out << "\\hline" << std::endl;
    out << "\\hline" << std::endl;
    out << "\\end{tabular}" << std::endl;
    out << "\\end{center}"  << std::endl;
    out << "\\end{table}"  << std::endl;   

 
    out.close();
    
    
    
    
    
    return;
}
void MiniTreeAnalyzer::PrintInvertedTableYields(){

 TString output=m_tables_dir;

  system("mkdir -p "+output);

  std::ofstream out;
  TString filename=output+"/Table_InvertedYields.tex";
  out.open(filename,std::ofstream::out);

    out << "\\begin{table}[htdp]"  << std::endl;
    out << "\\begin{center}"  << std::endl;

    out << "\\begin{tabular}{l";
    cout.flush();

       for (unsigned int process=0; process< PhysProcess.size(); process++){
		if (process==PhysProcess.size()-1) out << "|c}" << std::endl;
                else out << "|c";
	}

    out << "\\hline" << std::endl;
    out << "\\hline" << std::endl;
    out << " Cut " ;
    cout.flush();

       for (unsigned int process=0; process< PhysProcess.size(); process++){

		TString local=PhysProcess[process]->GetLatex().ReplaceAll("#","\\");
		if (local.Contains("\\")) local="$"+local+"$";

                if (process==PhysProcess.size()-1) out << " & " << local << " \\\\" << std::endl;
                else out << " & " << local;
       } 
 
    for (unsigned int fs=0; fs<channels.size(); fs++){


	out << "$"+channels[fs]->GetLatex().ReplaceAll("#","\\")+"$";
        cout.flush();

	for (unsigned int process=0; process< PhysProcess.size(); process++){

		float yield=PhysProcess[process]->GetYield(channels[fs]->GetName());
		float unc=PhysProcess[process]->GetStatistical(channels[fs]->GetName());

		if (process==PhysProcess.size()-1) out << " & $"<< std::fixed << std::setprecision(3) << yield << "\\pm"  << std::fixed << std::setprecision(3) << unc  << "$ \\\\" << std::endl;
		else out << " & $"<< std::fixed << std::setprecision(3) << yield << "\\pm"  << std::fixed << std::setprecision(3) << unc << "$";
	}

    } // end fs
  
    out << "\\hline" << std::endl;
    out << "\\hline" << std::endl;
    out << "\\end{tabular}" << std::endl;
    out << "\\end{center}"  << std::endl;
    out << "\\end{table}"  << std::endl;

    out.close();

 return;
}
void MiniTreeAnalyzer::PrintCutFlow(){

    std::cout << "#####################  Printing CUTFLOW  #################" << std::endl;    


    for (unsigned int fs=0; fs<channels.size(); fs++){

        std::cout << "- CutFlow for final state: " << channels[fs]->GetName() << std::endl;

        for (unsigned int process=0; process< PhysProcess.size(); process++){
    
              std::cout<< "-" << PhysProcess[process]->GetTitle() << ": " << std::endl;

              if (doPlotDistributions){
                 TCanvas *pC = new TCanvas;
                 pC->cd();
                 PhysProcess[process]->GetCutFlow(fs)->Draw("hist");
              }

              for (unsigned int cut=0; cut<PhysProcess[process]->GetCutFlow(fs)->GetNbinsX(); cut++){

                   std::cout << "--  " << PhysProcess[process]->GetCutFlow(fs)->GetXaxis()->GetBinLabel(cut+1) << " : "  << PhysProcess[process]->GetCutFlow(fs)->GetBinContent(cut+1) << std::endl;   

	      }


        }// process

    }//end fs

    
    std::cout << "##### WARNING : ADD TABLE IN LATEX WITH CUTFLOW WHEN YOU HAVE TIME!!!" << std::endl;


    
    return;
}
void MiniTreeAnalyzer::AddProcess(TString name,TString title,TString latex,int color,int line,TString process,bool isDD=false,TString cut="",TString weight="",TString tree="",int dict=0,TString syst_shape="",TString syst_norm=""){

    
    if (process != "isSig" && process != "isBkg" && process != "isData"){
    
        std::cout <<"- The following process '" << name << "' does not correspond to any of the configurations: 'isSig', 'isData' or 'isBkg'" << std::endl;
        std::cout <<"-- Please enter the right TString... ABORTING" << std::endl;
        exit(1);
    }
    
    if (process=="isSig") doSig=true;
    else if (process=="isData") doData=true;
    else if (process=="isBkg") {
       doBkg=true;
       //std::cout << "--- Setting doBkg to true " << std::endl;
    }
        
   std::vector<TString> samples=FillVector(name);
    
   
   if (process !="isSig"){ 
    PhysicsProcess* p_Temp  = new PhysicsProcess();
    
     for (unsigned int i=0; i<samples.size(); i++) {
      if (process=="isBkg") p_Temp->AddSample(samples[i],bkgDir,1);
      else if (process=="isSig") p_Temp->AddSample(samples[i],sigDir,1);
      else if (process=="isData") p_Temp->AddSample(samples[i],datDir,1);
     }

    
    p_Temp->SetStyle(name,title,latex,color,line,process); //kRed
    
    p_Temp->SetWeight(weight);
    p_Temp->SetTree(tree);
    p_Temp->SetCut(cut);
    p_Temp->SetDictionary(dict);   
    p_Temp->SetIsDataDriven(isDD);    
    p_Temp->SetSystematicVariations(syst_shape,syst_norm);


    PhysProcess.push_back(p_Temp);
   }
   
   
   if (process=="isSig"){
       for (unsigned int i=0; i<samples.size(); i++) {
    
	    std::string signal_input(samples[i].Data());	
	    std::vector<std::string> gettokens;	


	    tokenizeString(signal_input,' ',gettokens);


	    TString input_file;
	    input_file=gettokens[0];
	    TString local_tree=tree;	
	    std::vector<TString> grid;


	    if (gettokens.size() > 1) {
                local_tree=gettokens[1]; 	  
	         //test
         	//grid=GetMassPoints(local_tree,3,4);
            }
            //test
	    //else grid=GetMassPoints(input_file,6,7);
	    


            PhysicsProcess* p_Temp_sig  = new PhysicsProcess();
	    p_Temp_sig->AddSample(input_file,sigDir,1);
	    

            //test
	    TString name="Signal"; //"Signal_"+grid[0]+"_"+grid[1];
            if (gettokens.size() > 1) name=local_tree;

            //test 
	    TString latex_local=latex;//"mH"; //"m(#tilde#chi_{1}^{#pm}/#tilde#chi_{2}^{0},#tilde#chi_{1}^{0})=("+grid[0]+","+grid[1]+") GeV";
            //TString latex="("+grid[0]+","+grid[1]+")";	    

	    p_Temp_sig->SetStyle(name,name,latex_local,color-i,line,process);
	    p_Temp_sig->SetWeight(weight);
	    p_Temp_sig->SetCut(cut);
            p_Temp_sig->SetTree(local_tree);   
            p_Temp_sig->AddPositionsMass(grid);
	    p_Temp_sig->SetDictionary(dict);
            p_Temp_sig->SetIsDataDriven(isDD);
	    p_Temp_sig->SetSystematicVariations(syst_shape,syst_norm);

            PhysProcess.push_back(p_Temp_sig);
    
        
	}//end of samples
    
   }//adding signal samples
    
    
    return;
}

std::vector<TString> MiniTreeAnalyzer::GetMassPoints(TString filename,int x,int y){

  std::string sample(filename.Data());
  
  std::vector<std::string> gettokens;
  
  tokenizeString(sample, '_', gettokens); // only 6 and 7 components matter
  
  
  std::string mass1=gettokens[x];
  std::string mass2=gettokens[y];
  
  
  mass1.erase(mass1.size()-2,2);
  mass2.erase(mass2.size()-2,2);
  
  
  TString m1(mass1);
  TString m2(mass2);
  
  
  std::vector<TString> mass;
  
  mass.push_back(m1);
  mass.push_back(m2);
  
  
  return mass;
}
void MiniTreeAnalyzer::SettingConstants(){

    if (channels.size()==0){
    
        Channel *pCh= new Channel("All","All","");
        channels.push_back(pCh);
    
    }
    
    for (unsigned int i=0; i<channels.size(); i++){
     
        final_state.push_back(channels[i]->GetName());
        fs_latex.push_back(channels[i]->GetLatex());
        
    }
    
    
    return;
 

}

void MiniTreeAnalyzer::SaveHistos(bool save,TString folder="Histos"){

  m_save=save;
 
  m_histos_dir=folder;

  return;
}
void MiniTreeAnalyzer::SaveYieldsTables(bool save,TString folder="Tables"){

  m_save_tables=save;
  m_tables_dir=folder;


 return;
}
void MiniTreeAnalyzer::SavePlots(bool save,TString folder="Plots"){

 m_save_plots=save;
 m_plots_dir=folder;


 return;
}
void MiniTreeAnalyzer::AddVariable(VariableDistr* var){
    
    // std::cout << "-- Adding Variable .." << std::endl;
    
    variable.push_back(var);
    
    return;
}
void MiniTreeAnalyzer::AddVariable2D(VariableDistr* var_x,VariableDistr* var_y){

   VariableDistr2D *var = new VariableDistr2D(var_x,var_y);

   variable2D.push_back(var);
 
  return;
}




void MiniTreeAnalyzer::LoopOverSamples()
{
    if (m_debug>0) std::cout<<"##### In LoopOverSamples()" << std::endl;
   
    for (unsigned int i=0; i<PhysProcess.size(); i++){
        
        vector< vector<TH1F*> > output_histos_process;
        output_histos_process.clear();
        
        std::map<TString,vector< vector<TH1F*> > > output_histos_process_systematic;


	vector< vector<TH2F*> > output_histos_process_2D;
        output_histos_process_2D.clear();

	
        std::vector<TH1F*> output_histos_process_cutflow;
        output_histos_process_cutflow.clear();
        
	
        
        if (m_debug>0) cout << "#### STARTING PROCESS:"<< PhysProcess[i]->GetTitle() << endl;
        
        for (unsigned int j=0; j<PhysProcess[i]->Samples.size(); j++){
            
            //  cout << " -- background sample number : " << j << endl;
            vector< vector<TH1F*> > output_histos_sample;
            output_histos_sample.clear();
	    
	    vector< vector<TH2F*> > output_histos_sample_2D;
            output_histos_sample_2D.clear();

            
            vector<TH1F*> output_histos_sample_cutflow;
            output_histos_sample_cutflow.clear();
 
            std::map<TString,vector< vector<TH1F*> > > output_histos_sample_systematic;
           
            
            output_histos_sample=ReadInputFile(PhysProcess[i]->Samples[j]->GetFile(),PhysProcess[i]->GetTitle(),
                                               PhysProcess[i]->Samples[j]->GetScaleFactor(),
                                               PhysProcess[i]->isDataSigBkg(),PhysProcess[i]->GetIsDataDriven(),
					       PhysProcess[i]->GetTree(),PhysProcess[i]->GetWeight(),PhysProcess[i]->GetCut(),
                                               PhysProcess[i]->GetDictionary(),output_histos_sample_cutflow);
                
            
	    
	    if (variable2D.size()>0) output_histos_sample_2D=ReadInputFile2D(PhysProcess[i]->Samples[j]->GetFile(),PhysProcess[i]->GetTitle(),
                                            		  PhysProcess[i]->Samples[j]->GetScaleFactor(),
                                            		  PhysProcess[i]->isDataSigBkg(),PhysProcess[i]->GetIsDataDriven(),PhysProcess[i]->GetTree(),PhysProcess[i]->GetWeight(),PhysProcess[i]->GetCut(),PhysProcess[i]->GetDictionary());


	    std::map<TString,TString> systematics_list=PhysProcess[i]->GetSystematicsVariations();              

            if (systematics_list.empty() && !PhysProcess[i]->GetIsDataDriven() && PhysProcess[i]->isDataSigBkg()!="isData") systematics_list=GetSystematicsVariations();



            if (doSystematicVariations){

  		   
		   for (std::map<TString,TString>::iterator it=systematics_list.begin(); it!=systematics_list.end(); ++it){

			TString local_tree=PhysProcess[i]->GetTree();
			TString local_weight=PhysProcess[i]->GetWeight();

			if (it->second=="isShape") local_tree.ReplaceAll("nom",it->first); //PhysProcess[i]->GetTitle()+"_"+it->first;
		        else if (it->second=="isNorm") local_weight=it->first;

			output_histos_sample_systematic[it->first] = ReadInputFile(PhysProcess[i]->Samples[j]->GetFile(),PhysProcess[i]->GetTitle(),
                                                                      PhysProcess[i]->Samples[j]->GetScaleFactor(),
                                               			      PhysProcess[i]->isDataSigBkg(),PhysProcess[i]->GetIsDataDriven(),
                                               			      local_tree,local_weight,PhysProcess[i]->GetCut(),
                                               			      PhysProcess[i]->GetDictionary(),output_histos_sample_cutflow);           

		   } //end for syst 

            } // end if doSyst

	    
            
            if (j==0) {
                output_histos_process=output_histos_sample;
                
                output_histos_process_cutflow=output_histos_sample_cutflow;
       
                 if (variable2D.size()>0) output_histos_process_2D=output_histos_sample_2D; 
	 
	         /*
                 for (unsigned int fs=0; fs<output_histos_process_cutflow.size(); fs++){
                  for (unsigned int var=0; var< output_histos_process[fs].size(); var++){
                        cout << "-- Title X-axis: " << output_histos_process[fs][var]->GetXaxis()->GetTitle() << endl;
                  }
                   std::cout << "--- The histogram has n bins: " << output_histos_sample_cutflow[fs]->GetNbinsX()  << std::endl;       
	         }*/
            }
            else {
             
	       for (unsigned int ff=0; ff<output_histos_sample.size(); ff++){
                    for (unsigned int var=0; var< output_histos_sample[ff].size(); var++){
                        //   cout << "-- Integral of final state: " << ff << ", var : " << var << " --> " << output_histos_sample[ff][var]->Integral() << endl;
                        //       if (j==0) output_histos_process[fs][var]=(TH1F*)output_histos_sample[fs][var]->Clone();
                        output_histos_process[ff][var]->Add(output_histos_sample[ff][var]);
                    } //end for var
                    
		    
		    if (variable2D.size()>0){
                    	for (unsigned int var_2D=0; var_2D< output_histos_sample_2D[ff].size(); var_2D++){
				output_histos_process_2D[ff][var_2D]->Add(output_histos_sample_2D[ff][var_2D]);
          	   	 }//end for var_2D
		     }

                    output_histos_process_cutflow[ff]->Add(output_histos_sample_cutflow[ff]);

		    
                } //end for ff

            }// end else
	    
	    
	    
	     if (doSystematicVariations){
		 
		   for (std::map<TString,TString>::iterator it=systematics_list.begin(); it!=systematics_list.end(); ++it){

			
			if (j==0) output_histos_process_systematic[it->first]=output_histos_sample_systematic[it->first];
                        else {
				for (unsigned int ff=0; ff<output_histos_sample_systematic[it->first].size(); ff++){
             			       for (unsigned int var=0; var< output_histos_sample_systematic[it->first][ff].size(); var++){
                                                              output_histos_process_systematic[it->first][ff][var]->Add(output_histos_sample_systematic[it->first][ff][var]);
                                        } //end for var
				 }//end of ff

			}//end add histos



		   } //end for syst per sample

            } // end if doSyst
	    
	    
	    
	    

        } // end loop of samples for each process
       /* 
        for (unsigned int ff=0; ff<output_histos_process.size(); ff++){
            for (unsigned int var=0; var< output_histos_process[ff].size(); var++){
               // cout << "-- Integral of final state: " << final_state[ff] << ", var : " << variable[var]->GetTitle() << " --> " << output_histos_process[ff][var]->Integral() << endl;
                //cout << "-- Title X-axis: " << output_histos_process[fs][var]->GetXaxis()->GetTitle() << endl;
            }
        }*/
        
 

        PhysProcess[i]->AddSystematicVariations(output_histos_process_systematic);
	
	 

        PhysProcess[i]->AddHistos(output_histos_process);
	if (variable2D.size()>0) PhysProcess[i]->AddHistos2D(output_histos_process_2D);
        PhysProcess[i]->AddCutFlow(output_histos_process_cutflow);
        
    } // end process loop
 
    if (m_debug>0) std::cout << "##### Ending loop over samples " << std::endl;
    
    
    
   
    return;
}
void MiniTreeAnalyzer::SetSystematicVariations(bool doSyst,TString list_shapes,TString list_norm,TString theory,TString OutputDir){

   doSystematicVariations=doSyst;


 if (doSystematicVariations){
   std::vector<TString> line_shapes;
   std::vector<TString> line_norm;
 
   if (list_shapes != "") line_shapes=FillVector(list_shapes);
   if (list_norm != "" ) line_norm=FillVector(list_norm);

   std::vector<TString> shapes;
   std::vector<TString> norm;
  
   for (unsigned int i=0; i<line_shapes.size(); i++){
  
        std::vector<std::string> temp;
        std::string sample(line_shapes[i].Data());

        tokenizeString(sample,'=',temp);

        shapes.push_back(temp[0]);

        List_Systematics[temp[0]]=temp[1];   
   }

   for (unsigned int i=0; i<line_norm.size(); i++){

        std::vector<std::string> temp;
        std::string sample(line_norm[i].Data());

        tokenizeString(sample,'=',temp);

        norm.push_back(temp[0]);

        List_Systematics[temp[0]]=temp[1];
   }
   

   for (unsigned int i=0; i<shapes.size(); i++) m_Systematics[shapes[i]]="isShape";
   for (unsigned int j=0; j<norm.size(); j++) m_Systematics[norm[j]]="isNorm";
   
   OutputSystDir=OutputDir;
 }


  file_theory_uncertainties=theory;


  return;
}

std::vector<TString> MiniTreeAnalyzer::FillVector(TString name){

    std::vector<TString> samples;

    std::ifstream infile(name.Data());

    if (!infile){
        std::cout <<"-- File: '"<< name <<"' does not exist --- ABORTING " << std::endl;
        exit(1);
    }

    std::string line;


    while (std::getline(infile,line) ){

        if (line.empty()) continue;

        TString ss = line;

        samples.push_back(ss);
    }

 return samples;
}

vector< vector<TH2F*> > MiniTreeAnalyzer::ReadInputFile2D(TString pFile,TString samplename,float sf,TString isdata,bool isDD,TString tree,TString process_weight,TString process_cut,int dict){
 
 	
    TFile *pInput = new TFile(pFile);
    if (!pInput){
        printf("-- File %s is missing\n", pInput->GetName());
        exit(1);
    }
    else {
        if (m_debug>0) cout << "-- Opening file :" << pFile << endl;
    }	


    TString tree_process=tree;
    if (tree=="") tree_process=m_tree;


    TTree* local_Tree = (TTree*)pInput->Get(tree_process);
    
    if (m_debug>0) std::cout<< "-- Tree entries:" << local_Tree->GetEntries() << std::endl;
    
    TH1F *pHfile=(TH1F*)pInput->Get("hCutFlow");

    int bin_xs=-1; //pHfile->GetXaxis()->FindBin("xsec");
    int bin_sumW=-1; //pHfile->GetXaxis()->FindBin("aSumW");
   
    

    float xs=1; //pHfile->GetBinContent(bin_xs);
    float sumW=1; //pHfile->GetBinContent(bin_sumW);


    float lumi=1;

   /*
    if (isdata != "isData" && !isDD && dict ==0){
        lumi=m_luminosity;
        bin_xs=pHfile->GetXaxis()->FindBin("xsec");
        bin_sumW=pHfile->GetXaxis()->FindBin("aSumW");
        xs=pHfile->GetBinContent(bin_xs);
        sumW=pHfile->GetBinContent(bin_sumW);
    }
   */

    if (xs <=0 || sumW <=0){
     std::cout << "--- This sample does not have cross-section or sumW --> Aborting!" << std::endl;
     std::cout << "---- sample: " << pFile << std::endl;
     std::cout << "---- xs = " << xs << std::endl;
     std::cout << "---- sumW = " << sumW << std::endl;
   
     //if (xs <=0 ) xs=1;
      exit(1);
    }

    float xsec_scale=xs*lumi/sumW;


    vector< vector<TH2F*> > histos_input_file;

    histos_input_file.clear();
	

    for (unsigned int ff=0; ff<final_state.size(); ff++){
        
        if (m_debug>0) std::cout<<"### In final state: " << final_state[ff] << std::endl;


        std::vector<EventCut*> local_cuts;
        local_cuts.clear();      
  

	if (channels[ff]->GetSampleCut(samplename) !="") {
	        
                 EventCut *pTempCut1= new EventCut(channels[ff]->GetSampleCut(samplename),samplename+": Sorry! NO_CUTFLOW_HERE");

                 local_cuts.push_back(pTempCut1);
        }
	else if (process_cut!=""){
	         
		 EventCut *pTempCut2 = new EventCut(process_cut,samplename+": Sorry! NO_CUTFLOW_HERE");
		 
                 local_cuts.push_back(pTempCut2);
        }
        else {

        	if (GetPreselectionCut()!=0){

	       		EventCut *pPreselection= new EventCut(GetPreselectionCut()->GetCut(),GetPreselectionCut()->GetName());
               		local_cuts.push_back(pPreselection);
       		}
            
            
       		if (channels[ff]->GetChannelCut()->GetCut()!=""){
            

            		EventCut *pChannelCut= new EventCut(channels[ff]->GetChannelCut()->GetCut(),channels[ff]->GetChannelCut()->GetName());

            		local_cuts.push_back(pChannelCut);
         
            		if (channels[ff]->GetCut().size()>0){

           			for (unsigned int cc=0; cc<channels[ff]->GetCut().size(); cc++){

           	     			EventCut *pChannelExtra= new EventCut(channels[ff]->GetCut()[cc]->GetCut(),channels[ff]->GetCut()[cc]->GetName());

               	     			local_cuts.push_back(pChannelExtra);
               			}
           		 }
            
       		 }//adding cuts for channel



        
        	for (unsigned int icuts=0; icuts<v_cuts.size(); icuts++){
 
 	    		EventCut *pTempCut= new EventCut(v_cuts[icuts]->GetCut(),v_cuts[icuts]->GetName());	
            		local_cuts.push_back(pTempCut);
        	} // /adding extra cuts
         		
 
        } // end adding cuts
        

	
        
       

        TString temp_cut="";

        for (unsigned int icuts=0; icuts<local_cuts.size(); icuts++){
                           
            if (icuts==0) temp_cut=local_cuts[icuts]->GetCut();
            else temp_cut=temp_cut+" && "+local_cuts[icuts]->GetCut();
            
            local_cuts[icuts]->SetCut(temp_cut);
            
            //std::cout << "-- Adding cuts: " << temp_cut << std::endl;
        }


        TString FullCut="";

        if (local_cuts.size()>0) FullCut=local_cuts[local_cuts.size()-1]->GetCut();
                
        
        TString weight=GetWeight();
        if (process_weight!="") weight=process_weight;
        if (channels[ff]->GetSampleWeight(samplename) !="") weight=channels[ff]->GetSampleWeight(samplename);

        if (isdata == "isData") weight="";

        
    
        if (weight !="" && FullCut != "" ) FullCut=weight+"*("+FullCut+")";
        else if (weight !="") FullCut=weight;
        
        if (m_debug>0) std::cout << "-- Applying weight:"<< weight << std::endl;
        if (m_debug>0) std::cout << "-- Full cut:" << FullCut  << std::endl;
        
        
        
        vector<TH2F*> pHistoVar;
	
	
	for (unsigned int i=0; i<variable2D.size(); i++){


             TH2F *pH = CreateHisto2D(variable2D[i]->GetVarX(),variable2D[i]->GetVarY(),final_state[ff]);
             //pH->Sumw2();

             TString Local_Variable_X=variable2D[i]->GetVarX()->GetName();
             TString Local_Variable_Y=variable2D[i]->GetVarY()->GetName(); 

             if (dict !=0) {
                Local_Variable_X=ConvertVariable(Local_Variable_X,dict);
                Local_Variable_Y=ConvertVariable(Local_Variable_Y,dict);
             }


	     local_Tree->Draw(Local_Variable_Y+":"+Local_Variable_X+" >> "+ pH->GetName(),FullCut,"goff");	


	     if (isdata!="isData") pH->Scale(xsec_scale);
             pH->SetDirectory(0);	

	     if (m_debug>0) std::cout << "-- ReadOpenFile2D(): " << variable2D[i]->GetVarX()->GetName() << " vs " << variable2D[i]->GetVarY()->GetName()  <<",  Integral -> " << pH->Integral() << std::endl;



             //if (printLog) std::cout << "___ Events will be normalized..." << std::endl;
             //pH->Scale(100.0/pH->Integral());

	     pHistoVar.push_back(pH);
                    
  
        }// end of var
        
        
        
        histos_input_file.push_back(pHistoVar);
 

        pHistoVar.clear();

        local_cuts.clear();
        local_cuts.shrink_to_fit();
        
    } // end for fs

    

    pInput->Close();
    


    return histos_input_file;


}
vector < vector<TH1F*> > MiniTreeAnalyzer::ReadInputFile(TString pFile,TString samplename,float sf,TString isdata,bool isDD,TString tree,TString process_weight,TString process_cut,int dict,std::vector<TH1F*>& cutflow){
    // Each PhysicsProcess is looping over a variable and over a final state.
    // TH1::AddDirectory(kFALSE);
   
    
    TString tree_process=tree;
    if (tree=="") tree_process=m_tree;
  
    
    TFile *pInput = new TFile(pFile);
    if (!pInput){
        printf("-- File %s is missing\n", pInput->GetName());
        exit(1);
    }
    else {

       if (m_debug>0) cout << "-- Opening file :" << pFile << " for "<< samplename << ", tree_name: "<< tree_process <<  endl;
    
    }



    TTree* local_Tree = (TTree*)pInput->Get(tree_process);
    
    if (m_debug>0) std::cout<< "-- Tree entries:" << local_Tree->GetEntries() << std::endl;
    
    TH1F *pHfile=(TH1F*)pInput->Get("hCutFlow");

 
    int bin_xs=-1; //pHfile->GetXaxis()->FindBin("xsec");
    int bin_sumW=-1; //pHfile->GetXaxis()->FindBin("aSumW");


    float xs=1; //pHfile->GetBinContent(bin_xs);
    float sumW=1; //pHfile->GetBinContent(bin_sumW);
    float lumi=1;

  /*
    if (isdata != "isData" && !isDD && dict ==0 ){
        lumi=m_luminosity;
        bin_xs=pHfile->GetXaxis()->FindBin("xsec");
        bin_sumW=pHfile->GetXaxis()->FindBin("aSumW");
        xs=pHfile->GetBinContent(bin_xs);
        sumW=pHfile->GetBinContent(bin_sumW);
    }
  */

    if (xs <=0 || sumW <=0){
     std::cout << "--- This sample does not have cross-section or sumW --> Aborting!" << std::endl;
     std::cout << "---- sample: " << pFile << std::endl;
     std::cout << "---- xs = " << xs << std::endl;
     std::cout << "---- sumW = " << sumW << std::endl;

     //if (xs <=0 ) xs=1;
          exit(1);
     }
     
    float xsec_scale=xs*lumi/sumW;

    /*std::cout << "--- bin content sumW : " << pHfile->GetBinContent(2) << std::endl;
    std::cout << "---- xs = " << xs << std::endl;
    std::cout << "---- sumW = " << sumW << std::endl;
    std::cout << "---- lumi  = " << m_luminosity << std::endl;
    std::cout << "---- Scaling by  : " << xsec_scale << std::endl;  
    */    

    vector< vector<TH1F*> > histos_input_file;
    
    histos_input_file.clear();
  
   
   
    for (unsigned int ff=0; ff<final_state.size(); ff++){
     
      //  std::cout<<"### In final state: " << final_state[ff] << std::endl;
   
        std::vector<EventCut*> local_cuts;
        local_cuts.clear();      
         

        if (channels[ff]->GetSampleCut(samplename) !="") {
	        
                 EventCut *pTempCut1= new EventCut(channels[ff]->GetSampleCut(samplename),samplename+": Sorry! NO_CUTFLOW_HERE");

                 local_cuts.push_back(pTempCut1);
        }
	else if (process_cut!=""){
	         
		 EventCut *pTempCut2 = new EventCut(process_cut,samplename+": Sorry! NO_CUTFLOW_HERE");
		 
                 local_cuts.push_back(pTempCut2);
        }
        else {

        	if (GetPreselectionCut()!=0){

	       		EventCut *pPreselection= new EventCut(GetPreselectionCut()->GetCut(),GetPreselectionCut()->GetName());
               		local_cuts.push_back(pPreselection);
       		}
            
            
       		if (channels[ff]->GetChannelCut()->GetCut()!=""){
            

            		EventCut *pChannelCut= new EventCut(channels[ff]->GetChannelCut()->GetCut(),channels[ff]->GetChannelCut()->GetName());

            		local_cuts.push_back(pChannelCut);
         
            		if (channels[ff]->GetCut().size()>0){

           			for (unsigned int cc=0; cc<channels[ff]->GetCut().size(); cc++){

           	     			EventCut *pChannelExtra= new EventCut(channels[ff]->GetCut()[cc]->GetCut(),channels[ff]->GetCut()[cc]->GetName());

               	     			local_cuts.push_back(pChannelExtra);
               			}
           		 }
            
       		 }//adding cuts for channel



        
        	for (unsigned int icuts=0; icuts<v_cuts.size(); icuts++){
 
 	    		EventCut *pTempCut= new EventCut(v_cuts[icuts]->GetCut(),v_cuts[icuts]->GetName());	
            		local_cuts.push_back(pTempCut);
        	} // /adding extra cuts
         		
 
        } // end adding cuts
        

        TString weight=GetWeight();
        if (process_weight!="") weight=process_weight;
        if (channels[ff]->GetSampleWeight(samplename) !="") weight=channels[ff]->GetSampleWeight(samplename);

        if (isdata == "isData") weight="";

        EventCut *pInput = new EventCut("","Input");
	local_cuts.insert(local_cuts.begin(),pInput);


        
        TH1F *pLocal = new TH1F("h_cutflow_"+final_state[ff],"h_cutflow_"+final_state[ff],local_cuts.size(),0,local_cuts.size());
        pLocal->Sumw2();
        pLocal->GetXaxis()->SetTitle("Cut");
        pLocal->GetYaxis()->SetTitle("Events");        


        TString temp_cut="";

        for (unsigned int icuts=0; icuts<local_cuts.size(); icuts++){
     
            pLocal->GetXaxis()->SetBinLabel(icuts+1,local_cuts[icuts]->GetName());
            
            if (icuts==0 || temp_cut=="") temp_cut=local_cuts[icuts]->GetCut();
            else temp_cut=temp_cut+" && "+local_cuts[icuts]->GetCut();
            
            local_cuts[icuts]->SetCut(temp_cut);
            
            //std::cout << "-- Adding cuts: " << temp_cut << std::endl;
        }


        TString FullCut="";

        if (local_cuts.size()>0) FullCut=local_cuts[local_cuts.size()-1]->GetCut();
                
    
        if (weight !="" && FullCut != "" ) FullCut=weight+"*("+FullCut+")";
        else if (weight !="") FullCut=weight;
        
        if (m_debug>0) std::cout << "-- Applying weight:"<< weight << std::endl;
        if (m_debug>0) std::cout << "-- Full cut:" << FullCut  << std::endl;
        
        
        vector<TH1F*> pHistoVar;

        
        for (unsigned int i=0; i<variable.size(); i++){
            
            TH1F *pH = CreateHisto(variable[i],final_state[ff]);
 

            TString Local_Variable=variable[i]->GetName();
            if (dict !=0) Local_Variable=ConvertVariable(variable[i]->GetName(),dict); 
          
            if (variable[i]->GetYield() && doCutflow){
                

                for (unsigned int cut_i=0; cut_i<local_cuts.size(); cut_i++){
                
                    TString temp_cutflow=local_cuts[cut_i]->GetCut();
 
                    if (weight !="" && temp_cutflow !="") temp_cutflow=weight+"*("+temp_cutflow+")";
                    else if (weight !="") temp_cutflow=weight;                   

                    local_Tree->Draw(Local_Variable+" >> "+ pH->GetName(),temp_cutflow,"goff");
 
                    if (isdata != "isData") pLocal->SetBinContent(cut_i+1,pH->Integral()*xsec_scale);
                    else pLocal->SetBinContent(cut_i+1,pH->Integral());
        
                    pLocal->SetDirectory(0);

	         } // end fill bin cutflow

            }//end do cutflow
            

            local_Tree->Draw(Local_Variable+" >> "+ pH->GetName(),FullCut,"goff");
            
            
            //if (isdata!="isData") scale_factor=scale_factor*lumi_correction;
            
            if (isdata!="isData") pH->Scale(xsec_scale);
            pH->SetDirectory(0);
            //float low=pH->GetBinLowEdge(1);
            //float up=pH->GetBinLowEdge(pH->GetNbinsX()+1);
            
            
            if (m_debug>0) std::cout << "-- ReadOpenFile(): " << variable[i]->GetTitle() << ",  Integral -> " << pH->Integral() << std::endl;

          //std::cout << "-- ReadOpenFile(): " << pH->GetBinLowEdge(1) << ", " << pH->GetBinLowEdge(pH->GetNbinsX()+1) << " :  " << pH->Integral(1,pH->GetNbinsX()) << std::endl;
           
           
            
            
            
           pHistoVar.push_back(pH);
           
         
  
        }// end of var
        
        
        
        histos_input_file.push_back(pHistoVar);
 

        if (doCutflow) cutflow.push_back(pLocal);        


        pHistoVar.clear();

        local_cuts.clear();
        local_cuts.shrink_to_fit();
        
    } // end for fs

    


    pInput->Close();
    



    return histos_input_file;
}
TString MiniTreeAnalyzer::ConvertVariable(TString var, int dict){

 TString new_var=var;

 TString toGeV="*1";

if (dict==1){

 if (var=="@jets.size()") new_var="Njet20";
 else if (var=="fabs(leps[0].eta-leps[1].eta)") new_var="fabs(DEtall)";
 //else if (var=="sig.MetRel") new_var="RelMET"+toGeV;
 else if (var=="sig.HT+sig.Met") new_var="Meff"+toGeV;
 else if (var=="sig.mlj")  new_var="Mlj"+toGeV;
 else if (var=="sig.Met") new_var="MET"+toGeV;
 else if (var=="sig.mT2") new_var="Mt2"+toGeV;
 else if (var=="leps[0].pt") new_var="lep_1stPt"+toGeV;
 else if (var=="leps[1].pt") new_var="lep_2ndPt"+toGeV;
 //else if (var.Contains("?")) new_var="Maxmt"+toGeV;
 else if (var=="l12.pt") new_var="Ptll"+toGeV;
 else if (var=="l12.m") new_var="Mll"+toGeV;
 else if (var=="TMath::Sqrt(2*leps[0].pt*sig.Met*(1-TMath::Cos(leps[0].MET_dPhi)))") new_var="Mt"+toGeV;

}
else if (dict==2){

 toGeV="*0.001";

 if (var=="@jets.size()") new_var="nJets25";
 else if (var=="fabs(leps[0].eta-leps[1].eta)") new_var="fabs(DeltaEtaLep)";
 //else if (var=="sig.MetRel") new_var="RelMET"+toGeV;
 else if (var=="sig.HT+sig.Met") new_var="meff"+toGeV;
 else if (var=="sig.mlj")  new_var="mljj_comb"+toGeV;
 else if (var=="sig.Met") new_var="met"+toGeV;
 else if (var=="sig.mT2") new_var="MT2"+toGeV;
 else if (var=="leps[0].pt") new_var="Pt_l"+toGeV;
 else if (var=="leps[1].pt") new_var="Pt_subl"+toGeV;
 //else if (var.Contains("?")) new_var="Maxmt"+toGeV;
 //else if (var=="l12.pt") new_var="Ptll"+toGeV;
 //else if (var=="l12.m") new_var="Mll"+toGeV;
 else if (var=="TMath::Sqrt(2*leps[0].pt*sig.Met*(1-TMath::Cos(leps[0].MET_dPhi)))") new_var="mt"+toGeV;


}

else if (dict==3){
  toGeV="*0.001";

  if (var=="met") new_var="met"+toGeV;
  else if (var=="meff") new_var="meff"+toGeV;
  

}

 return new_var;
}
void MiniTreeAnalyzer::PrintStatisticalUncertainty(){

 for (unsigned int var=0; var< variable.size(); var++){
        if (variable[var]->GetYield()==true){

            for (unsigned int process=0; process< PhysProcess.size(); process++){

                 if (PhysProcess[process]->isDataSigBkg()=="isData") continue;
                 for (unsigned int ff=0; ff<final_state.size(); ff++){

                     TH1F *pH=(TH1F*)PhysProcess[process]->GetHistos(ff,var)->Clone();

                     double uncertainty=0;
                     double yield=0;

                     yield=pH->IntegralAndError(1,pH->GetNbinsX(),uncertainty);                   

                     //std::cout << PhysProcess[process]->GetTitle() << ", " << final_state[ff] << ", " << uncertainty << std::endl;
                    
                
                } // end for process

            } // end for fs
        } // end if var==yield
 } // end for var
         

 return;
}
void MiniTreeAnalyzer::ComputeYields(){
    
    
    for (unsigned int var=0; var< variable.size(); var++){
        if (variable[var]->GetYield()==true){
            
            for (unsigned int ff=0; ff<final_state.size(); ff++){
                
              
                THStack *p_bkg= new THStack("bkg_total"+final_state[ff]+"_"+variable[var]->GetName(),"");
                
                if (m_debug>-1) std::cout << "###########  Yields for final state: " << final_state[ff] << std::endl;
              
                 
                for (unsigned int process=0; process< PhysProcess.size(); process++){
               
                 
     
                    TH1F *pH=(TH1F*)PhysProcess[process]->GetHistos(ff,var)->Clone();
                    
                    double uncertainty=0;
                    double yield=0;
                    
                    yield=pH->IntegralAndError(1,pH->GetNbinsX(),uncertainty);
                    
                 
                       
                    Yields *p_yield= new Yields(yield,uncertainty);
                    
 		
		    TString local_file_theory_uncertainties=file_theory_uncertainties;

                    PhysProcess[process]->SetYield(final_state[ff],p_yield);
		    if (doSystematicVariations) PhysProcess[process]->SetSystematic(ff,var,final_state[ff],local_file_theory_uncertainties); 
		    
		     std::map<TString,vector< vector<TH1F*> >> LocalMap=PhysProcess[process]->GetMapSystematics();
		    

		    if (m_debug>-1) {
                         if (doSystematicVariations) std::cout << "ComputeYields()::" << PhysProcess[process]->GetTitle() << " ->  yield:" << yield << " +/- " << uncertainty << " +" << PhysProcess[process]->GetSystematic(final_state[ff])[0] << " - " << PhysProcess[process]->GetSystematic(final_state[ff])[1] <<std::endl;			 
                         else std::cout << "ComputeYields()::" << PhysProcess[process]->GetTitle() << " ->  yield:" << yield << " +/- " << uncertainty << std::endl;
                     }
                   
                    //std::vector<double> Syst=PhysProcess[process]->GetTotalSystematic(final_state[fs]);
                    
                    //std::cout << "-----------------" << PhysProcess[process]->GetName() << " ->  yield:" << yield << " +/- "
                    //<< uncertainty << " : -" << Syst[0] << " / +" << Syst[1] << std::endl;
                    
                      
                    if (PhysProcess[process]->isDataSigBkg()=="isBkg") p_bkg->Add((TH1F*)pH->Clone());
                    else if (PhysProcess[process]->isDataSigBkg()=="isData") TotalData[final_state[ff]]=p_yield;
                    else if (PhysProcess[process]->isDataSigBkg()=="isSig") TotalSig[PhysProcess[process]->GetName()][final_state[ff]]=p_yield;
                    
                } // end for process
                
               
                if (doBkg){
                    double total_bkg_uncertainty=0;
                    double total_bkg=0;
                    
                    TH1F *h_base=(TH1F*)p_bkg->GetStack()->Last();
                    
                    total_bkg=h_base->IntegralAndError(1,h_base->GetNbinsX(),total_bkg_uncertainty);
                    
                    Yields *p_bkg_total= new Yields(total_bkg,total_bkg_uncertainty);
                    
                    TotalBkg[final_state[ff]]=p_bkg_total;
       
		    if (doSystematicVariations) ComputeSystematicVariations();
             
                    if (m_debug>-1){
                        if (doSystematicVariations) std::cout << "--- Total background : "<< TotalBkg[final_state[ff]]->GetYield() << " +/- " << TotalBkg[final_state[ff]]->GetStatistical()  << " + " << TotalBkg[final_state[ff]]->GetSystematic()[0] << " - " << TotalBkg[final_state[ff]]->GetSystematic()[1] <<  std::endl;
                    else std::cout << "--- Total background : "<< TotalBkg[final_state[ff]]->GetYield() << " +/- " << TotalBkg[final_state[ff]]->GetStatistical() << std::endl;
                 
                     }
    /*
		    if (doSystematicVariations){

                       std::vector<float> syst=
		    }*/

                } // end if doBkg
                
             

            } // end for fs
        } // end if var==yield
    } // end for var
    
    
    
    if (doSignificance) GetSignificance();
    
    
    return;
}

TH1F* MiniTreeAnalyzer::CreateHisto(VariableDistr*  var,TString ch){
    
    TH1F *pH; // = new TH1F("hist"+var->GetTitle()+"_"+fs,";"+var->GetXLabel()+";"+var->GetYLabel(),var->GetNbins(),var->GetLowAxis(),var->GetUpAxis());

    TString ylabel=var->GetYLabel();

    if (m_weight=="") ylabel.ReplaceAll("Events","Unweighted events");
    if (NormalizeToUnit) ylabel="Normalized "+ylabel;

    
    if (var->IsArray()) {
        pH = new TH1F("hist"+var->GetTitle()+"_"+ch,"hist"+var->GetTitle()+"_"+ch,var->GetNbins(),var->GetBinning());
        
        pH->GetYaxis()->SetTitle(ylabel);
        pH->GetXaxis()->SetTitle(var->GetXLabel());
    }
    else pH = new TH1F("hist"+var->GetTitle()+"_"+ch,";"+var->GetXLabel()+";"+ylabel,var->GetNbins(),var->GetLowAxis(),var->GetUpAxis());
    
       
    pH->Reset();
    
    if (!NormalizeToUnit) pH->Sumw2();
    
    return pH;
}

TH2F* MiniTreeAnalyzer::CreateHisto2D(VariableDistr*  var_x,VariableDistr*  var_y, TString fs){

 TH2F *pH=0;


 if (var_x->IsArray() && !var_y->IsArray()){
     pH=new TH2F("hist"+var_x->GetTitle()+"_"+var_y->GetTitle()+"_"+fs,"hist"+var_x->GetTitle()+"_"+var_y->GetTitle()+"_"+fs,var_x->GetNbins(),var_x->GetBinning(),var_y->GetNbins(),var_y->GetLowAxis(),var_y->GetUpAxis());  
 }
 else if (!var_x->IsArray() && var_y->IsArray()){
     pH=new TH2F("hist"+var_x->GetTitle()+"_"+var_y->GetTitle()+"_"+fs,"hist"+var_x->GetTitle()+"_"+var_y->GetTitle()+"_"+fs, var_x->GetNbins(),var_x->GetLowAxis(),var_x->GetUpAxis(),var_y->GetNbins(),var_y->GetBinning());
 }
 else if (var_x->IsArray() && var_y->IsArray()){
     pH=new TH2F("hist"+var_x->GetTitle()+"_"+var_y->GetTitle()+"_"+fs,"hist"+var_x->GetTitle()+"_"+var_y->GetTitle()+"_"+fs,var_x->GetNbins(),var_x->GetBinning(),var_y->GetNbins(),var_y->GetBinning());
 }
 else {
    pH=new TH2F("hist"+var_x->GetTitle()+"_"+var_y->GetTitle()+"_"+fs,"hist"+var_x->GetTitle()+"_"+var_y->GetTitle()+"_"+fs,
		    var_x->GetNbins(),var_x->GetLowAxis(),var_x->GetUpAxis(),var_y->GetNbins(),var_y->GetLowAxis(),var_y->GetUpAxis());
 }

 pH->Reset();

 TString zlabel="Events";
 if (m_weight=="") zlabel.ReplaceAll("Events","Unweighted events");

 pH->GetZaxis()->SetTitle(zlabel);

 pH->GetXaxis()->SetTitle(var_x->GetXLabel());
 pH->GetYaxis()->SetTitle(var_y->GetXLabel());


 pH->Sumw2();

 return pH;
}


void MiniTreeAnalyzer::AddChannel(Channel *p_ch){
    
    channels.push_back(p_ch);
    
    return;
    
}
void MiniTreeAnalyzer::AddCut(TString cut,TString name){
    

    EventCut *pTemp= new EventCut(cut,name);
    
    v_cuts.push_back(pTemp);
    
    
    return;
}
void MiniTreeAnalyzer::AddPreselectionCut(TString cut,TString name){
    
    //std::cout << "-- Adding preselection cut :" << cut << std::endl;
    
    p_Cut_preselection= new EventCut(cut,name);
    
    return;
}
void MiniTreeAnalyzer::GetSignificance()
{
    if (!doSig){
        std::cout << "-- Signal samples have not been added -> Significance will not be compute! JUMPING " << std::endl;
        return;
    }
    
        
    for (auto it_sig: TotalSig){
       if (m_debug>0){
          std::cout << "--------------------------------------------------------------------------------" <<std::endl;
          std::cout << "####### Results for " << it_sig.first << " ############" << std::endl;
        }

         double sig2=0;


        for (auto it_ch: it_sig.second){
            
            double sig=it_ch.second->GetYield();
            double sig_unc=it_ch.second->GetStatistical();
            double bkg=TotalBkg[it_ch.first]->GetYield();
            double bkg_unc=TotalBkg[it_ch.first]->GetStatistical();
   
            double Root_sig=-1;
            double Null_sig=-1;
            double Stat_sig=-1;
            double Gen_sig=-1;
            
            // std::cout << it_sig.first << " : " << it_ch.first <<" : "<<  sig  << std::endl;
            // std::cout << "- Total background: " <<  bkg  << std::endl;


            if (bkg>0 && sig>0){
                 if (m_root_sig) Root_sig=RooStats::NumberCountingUtils::BinomialExpZ(sig,bkg,TMath::Sqrt(pow(bkg_unc/bkg,2)+pow(m_rel_syst,2)));       
                 if (m_standard) Null_sig=ComputeSignificance(sig,bkg);
                 if (m_stat_standard) Stat_sig=ComputeStatisticalSignificance(sig,bkg);
                 if (m_stat_forum) Gen_sig=ComputeGeneralSignificance(sig,bkg, TMath::Sqrt(pow(bkg_unc,2)+pow(m_rel_syst*bkg,2))); 
             }

     
            Significance[it_sig.first][it_ch.first]=Null_sig;
            StatisticalSignificance[it_sig.first][it_ch.first]=Stat_sig;
            GeneralSignificance[it_sig.first][it_ch.first]=Gen_sig;
            RootSignificance[it_sig.first][it_ch.first]=Root_sig;            
           
           if (m_debug>0) {
            std::cout << "------ Channel : " << it_ch.first << std::endl;
            std::cout << "Signal yield :" << sig << " +/- " << sig_unc << std::endl;
            std::cout << "Total bkg: " << bkg << " +/- " << bkg_unc << std::endl;
	    if (m_root_sig) std::cout << "-RootSignificance : " << Root_sig << std::endl;
            if (m_stat_forum) std::cout << "-Significance (statistical forum - recommended one!) : " << Gen_sig << std::endl;
            if (m_standard) std::cout << "-Significance vs. null hypothesis :" << Null_sig << std::endl;
            if (m_stat_standard) std::cout << "-Statatistical Significance: " << Stat_sig << std::endl;
            std::cout << "---" <<std::endl;
           }            

            sig2=sig2+pow(Root_sig,2);

        }

            double combined_sig=TMath::Sqrt(sig2);

	    if (m_debug>0) std::cout<< "--- Combined significance (with RootSignificance): " << combined_sig << std::endl;

            CombinedSignificance[it_sig.first]=combined_sig;


         //   std::cout << "--------------------------------------------------------------------------------" <<std::endl;
    }
    
    
    // for (auto it_sig:Significance){
    //  for (auto it_ch:it_sig.second){
    
    //       std::cout << it_sig.first << " : " << it_ch.first <<" : "<< it_ch.second << std::endl;     
    
    //  }
    // }
    
    
    return;
}
double MiniTreeAnalyzer::ComputeSignificance(double sig,double bkg)
{
    return sig/TMath::Sqrt(bkg);
}
double MiniTreeAnalyzer::ComputeStatisticalSignificance(double sig,double bkg)
{
    return sig/TMath::Sqrt(sig+bkg);
}
double MiniTreeAnalyzer::ComputeGeneralSignificance(double sig,double bkg,double ebkg)
{
   /*Significance recommended by Statistics Forum */

    double sigma2 = ebkg*ebkg;  
    double bkg2 = bkg*bkg;
    
    double f1 = sig*TMath::Log(sig*(bkg +sigma2)/(bkg2+sig*sigma2));
    double f2 = (bkg2/sigma2)*TMath::Log((bkg2+sig*sigma2)/(bkg*(bkg+sigma2)));
 
    double val=TMath::Sqrt(2*(f1-f2));

    if (sig < bkg) val=-val;

    return val;

    //return TMath::Sqrt(2*((sig+bkg)*TMath::Log(1+(sig/bkg))-sig));
}

std::vector<TString> MiniTreeAnalyzer::ParticleMCTC(TString var){

   
 std::vector<TString> sParticleType;

 if (var.Contains("truthType")){
   sParticleType.push_back("Unk.");
   sParticleType.push_back("Unk.El");
   sParticleType.push_back("Iso.El");
   sParticleType.push_back("NonIso.El");
   sParticleType.push_back("Bkg.El");
   sParticleType.push_back("Unk.Mu");
   sParticleType.push_back("Iso.Mu");
   sParticleType.push_back("NonIso.Mu");
   sParticleType.push_back("Bkg.Mu");
   sParticleType.push_back("Unk.Tau");
   sParticleType.push_back("Iso.Tau");
   sParticleType.push_back("NonIsoTau");
   sParticleType.push_back("BkgTau");
   sParticleType.push_back("UnknownPhoton");
   sParticleType.push_back("IsoPhoton");
   sParticleType.push_back("NonIsoPhoton");
   sParticleType.push_back("BkgPhoton");
   sParticleType.push_back("Hadron");
   sParticleType.push_back("Neutrino");
   sParticleType.push_back("NuclFrag");
   sParticleType.push_back("NonPrimary");
   sParticleType.push_back("GenParticle");
   sParticleType.push_back("SUSYParticle");
   sParticleType.push_back("BBbarMesonPart");
   sParticleType.push_back("BottomMesonPart");
   sParticleType.push_back("CCbarMesonPart");
   sParticleType.push_back("CharmedMesonPart");
   sParticleType.push_back("BottomBaryonPart");
   sParticleType.push_back("CharmedBaryonPart");
   sParticleType.push_back("StrangeBaryonPart");
   sParticleType.push_back("LightBaryonPart");
   sParticleType.push_back("StrangeMesonPart");
   sParticleType.push_back("LightMesonPart");
   sParticleType.push_back("BJet");
   sParticleType.push_back("CJet");
   sParticleType.push_back("LJet");
   sParticleType.push_back("GJet");
   sParticleType.push_back("TauJet");
   sParticleType.push_back("UnknownJet");
}
else if (var.Contains("truthOrig")){
   sParticleType.push_back("NonDef.");
   sParticleType.push_back("SingleElec");
   sParticleType.push_back("SingleMuon");
   sParticleType.push_back("SinglePhot");
   sParticleType.push_back("SingleTau");
   sParticleType.push_back("PhotonConv");
   sParticleType.push_back("DalitzDec");
   sParticleType.push_back("ElMagProc");
   sParticleType.push_back("Mu");
   sParticleType.push_back("TauLep");
   sParticleType.push_back("top");
   sParticleType.push_back("QuarkWeakDec");
   sParticleType.push_back("W");
   sParticleType.push_back("Z");
   sParticleType.push_back("Higgs");
   sParticleType.push_back("HiggsMSSM");
   sParticleType.push_back("WZMSSM");
   sParticleType.push_back("WBosonLRSM");
   sParticleType.push_back("NuREle");
   sParticleType.push_back("NuRMu ");
   sParticleType.push_back("NuRTau");
   sParticleType.push_back("LQ");
   sParticleType.push_back("SUSY");
   sParticleType.push_back("LightMeson");
   sParticleType.push_back("StrangeMeson");
   sParticleType.push_back("CharmedMeson");
   sParticleType.push_back("BottomMeson");
   sParticleType.push_back("CCbarMeson");
   sParticleType.push_back("JPsi");
   sParticleType.push_back("BBbarMeson");
   sParticleType.push_back("LightBaryon");
   sParticleType.push_back("StrangeBaryon");
   sParticleType.push_back("CharmedBaryon");
   sParticleType.push_back("BottomBaryon");
   sParticleType.push_back("PionDecay");
   sParticleType.push_back("KaonDecay");
   sParticleType.push_back("BremPhot");
   sParticleType.push_back("PromptPhot");
   sParticleType.push_back("UndrPhot");
   sParticleType.push_back("ISRPhot");
   sParticleType.push_back("FSRPhot");
   sParticleType.push_back("NucReact");
   sParticleType.push_back("PiZero");
   sParticleType.push_back("VV");
   sParticleType.push_back("ZorHeavyBoson");
   sParticleType.push_back("QCD");
}
  
  

  return sParticleType;
}
void MiniTreeAnalyzer::GetLabelMCTruthClassifier(TH2F* &pH,TString var_x,TString var_y){

 
 std::cout <<"####################################  Getting label for MCTruthClassifier -  2D" << std::endl;

  // Fill x-axis
 std::vector<TString> particle_x=ParticleMCTC(var_x);

 if (particle_x.size()!=pH->GetNbinsX()) return; //bins of histograms has to be compatible with vector size of Type or Origin.

 for (unsigned int i=0; i<pH->GetNbinsX(); i++){
  
   //std::cout<<"-- Integral x : " << i << "-- val: " << integral_x << std::endl; 

   //if (integral_x !=0 ){
       pH->GetXaxis()->SetBinLabel(i+1,particle_x[i]);
    //   std::cout << "---- label : " << particle_x[i] << std::endl;
  // }

 }


 if (var_x.Contains("truthType")) pH->GetXaxis()->SetRangeUser(1,7);

  
  // Fill y-axis
 std::vector<TString> particle_y=ParticleMCTC(var_y);

 if (particle_y.size()!=pH->GetNbinsY()) return; //bins of histograms has to be compatible with vector size of Type or Origin.

 for (unsigned int i=0; i<pH->GetNbinsY(); i++){
  
   //std::cout<<"-- Integral y : " << i << "-- val: " << integral_y << std::endl;

  // if (integral_y !=0 ){
       pH->GetYaxis()->SetBinLabel(i+1,particle_y[i]);
     //  std::cout << "---- label : " << particle_y[i] << std::endl;
  // }

 }

 //->GetXaxis()->SetRangeUser(1
 
 
 if (var_y.Contains("truthType")) pH->GetYaxis()->SetRangeUser(1,7);



std::cout<< "#########################################  Setting MCTRUTHCLASSIFIER LABELS:  DONE!!  #####" << std::endl;
 


 return;
}
void MiniTreeAnalyzer::GetLabelMCTruthClassifier(TH1F* &pH,TString var){

 std::cout<< "#########################################  Setting MCTRUTHCLASSIFIER LABELS #####" << std::endl;

 std::vector<TString> sParticleType=ParticleMCTC(var);


 if (sParticleType.size()!=pH->GetNbinsX()) return; //bins of histograms has to be compatible with vector size of Type or Origin.

 int bin=0;

 for (unsigned int i=0; i<pH->GetNbinsX(); i++){
   float content=pH->GetBinContent(i+1);
  
   if (content !=0 ){
       pH->GetXaxis()->SetBinLabel(i+1,sParticleType[i]);
       bin=i+1;
   }

 }

 //->GetXaxis()->SetRangeUser(1
 
 
 if (var.Contains("truthType")) pH->GetXaxis()->SetRangeUser(pH->GetXaxis()->GetBinLowEdge(1),pH->GetXaxis()->GetBinLowEdge(bin+1));



std::cout<< "#########################################  Setting MCTRUTHCLASSIFIER LABELS:  DONE!!  #####" << std::endl;


 return;
}
void MiniTreeAnalyzer::tokenizeString(std::string& str, char delim, std::vector<std::string>& tokens)
{

    //if (delim !=' ') str.erase(remove_if(str.begin(), str.end(), ::isspace), str.end()); 
    std::replace( str.begin(), str.end(), '\t', ' ');
    std::string str_local=regex_replace(str,regex("\\s{2,}")," ");

    tokens.clear();
    std::istringstream iss(str_local);
    std::string token;
    while ( std::getline(iss, token, delim) ){
        //boost::algorithm::trim(token);
        //token=token;
        //
        token.erase(remove_if(token.begin(), token.end(), ::isspace), token.end()); 
        //if (token==" ") continue;

        tokens.push_back(token);
    }
   return; 


}
void MiniTreeAnalyzer::ComputeSystematicVariations(){


    //get full vector of systematics:
    std::vector<TString> FullSystematics=GetFullListSystematics();

    std::map<TString,vector< vector<TH1F*> > > output_histos_total_per_systematic;
  


   //adding linearly same systematic for different samples since they are fully correlated
   //
   for (unsigned int syst=0; syst<FullSystematics.size(); syst++){

       for (unsigned int process=0; process< PhysProcess.size(); process++){


          if (PhysProcess[process]->isDataSigBkg()=="isBkg"){

               std::map<TString,vector< vector<TH1F*> > > syst_map=PhysProcess[process]->GetMapSystematics();

	      vector< vector<TH1F*> > local_syst;
               
               // check if systematic exist in sample, otherwise use nominal
              if (syst_map.count(FullSystematics[syst])>0) local_syst=syst_map[FullSystematics[syst]]; 
	      else local_syst=PhysProcess[process]->GetVectorHistos();


              vector< vector<TH1F*> > local_syst_clone;
	      //get clone of local_syst to avoid rewriting the histos. Subtracting nominal values to get directly the variation.
	       for (unsigned int ff=0; ff< final_state.size(); ff++){
                       
                                 vector<TH1F*> local_1;
			         for (unsigned int var=0; var<variable.size(); var++){         
      
					TH1F *pH_local=0;
					pH_local=(TH1F*)local_syst[ff][var]->Clone();
					pH_local->Add(PhysProcess[process]->GetHistos(ff,var),-1);

                                        local_1.push_back(pH_local);
          
                                 }// end ff
				 local_syst_clone.push_back(local_1);
                }//end var


              // check if systematic exist in total map, otherwise add variations linearly 
	      if (output_histos_total_per_systematic.count(FullSystematics[syst]) >0 ){

		       for (unsigned int var=0; var< variable.size(); var++){
                      		for (unsigned int ff=0; ff<final_state.size(); ff++){
                   
					output_histos_total_per_systematic[FullSystematics[syst]][ff][var]->Add(local_syst_clone[ff][var]);

                      	         }// end ff
                	}//end var
              } //end if syst exist total
	      else  output_histos_total_per_systematic[FullSystematics[syst]]=local_syst_clone;           

          }//end if bkg

        }// end process
     } // end syst


   //Split list in up and down
   std::map<TString,vector< vector<TH1F*> > > Syst_up;
   std::map<TString,vector< vector<TH1F*> > > Syst_down;

   for (unsigned int syst=0; syst<FullSystematics.size(); syst++){

        if (FullSystematics[syst].Contains("up") || FullSystematics[syst].Contains("Up") || FullSystematics[syst].Contains("UP") || FullSystematics[syst].Contains("1up")) Syst_up[FullSystematics[syst]]=output_histos_total_per_systematic[FullSystematics[syst]];
        else if (FullSystematics[syst].Contains("down") || FullSystematics[syst].Contains("Down") || FullSystematics[syst].Contains("DOWN") || FullSystematics[syst].Contains("1dn") ) Syst_down[FullSystematics[syst]]=output_histos_total_per_systematic[FullSystematics[syst]];
        else {
               Syst_up[FullSystematics[syst]]=output_histos_total_per_systematic[FullSystematics[syst]];
	       Syst_down[FullSystematics[syst]]=output_histos_total_per_systematic[FullSystematics[syst]];
        }

        if (FullSystematics[syst].Contains("JET_JER_SINGLE_NP")) Syst_down[FullSystematics[syst]]=output_histos_total_per_systematic[FullSystematics[syst]];

   }//end for list syst


   std::vector< std::vector<TH1F*> > syst_up2;
   std::vector< std::vector<TH1F*> > syst_down2; 



   //add quadrature of different syst.		
  int counter_up=0;
  for (std::map<TString, vector< vector<TH1F*> > >::iterator it=Syst_up.begin(); it!=Syst_up.end(); ++it){
 
 
        if (counter_up==0){
	     syst_up2=it->second;
	
   	     for (unsigned int ff=0; ff<final_state.size(); ff++){               
		 for (unsigned int var=0; var<variable.size(); var++){

                    syst_up2[ff][var]=GetSquareTH1F(it->second[ff][var]);
		 }
	     }

        }
        else {

           for (unsigned int ff=0; ff< final_state.size(); ff++){
   
                 for (unsigned int var=0; var<variable.size(); var++){


 		     syst_up2[ff][var]->Add(GetSquareTH1F(it->second[ff][var]));

                  } //end var
          }// end ff
    
       }

      counter_up=counter_up+1;
   } 

 
  int counter_down=0;
  for (std::map<TString, vector< vector<TH1F*> > >::iterator it=Syst_down.begin(); it!=Syst_down.end(); ++it){


        if (counter_down==0){
             syst_down2=it->second;

             for (unsigned int ff=0; ff<final_state.size(); ff++){

                 for (unsigned int var=0; var<variable.size(); var++){

                    syst_down2[ff][var]=GetSquareTH1F(it->second[ff][var]);
                 }
             }

        }
        else {

           for (unsigned int ff=0; ff< final_state.size(); ff++){

                 for (unsigned int var=0; var<variable.size(); var++){


                     syst_down2[ff][var]->Add(GetSquareTH1F(it->second[ff][var]));

                  } //end var
          }// end ff


         }

      counter_down=counter_down+1;
   }




 


 std::vector< std::vector<TH1F*> > total_syst_up;
 std::vector< std::vector<TH1F*> > total_syst_down;
 
  
 //Getting Sqrt(bin) to get final total systematic for bkg
  for (unsigned int ff=0; ff< final_state.size(); ff++){
 
         std::vector<TH1F*> local_down;
         std::vector<TH1F*> local_up;
         for (unsigned int var=0; var<variable.size(); var++){   

                TH1F *local_exp_up=(TH1F*)syst_up2[ff][var]->Clone();
                TH1F *local_exp_down=(TH1F*)syst_down2[ff][var]->Clone();
 
                TFile *pFile = new TFile("Shapes/TH1F_TotalUp_"+variable[var]->GetTitle()+"_"+final_state[ff]+".root");                


                TH1F *theory_up=(TH1F*)pFile->Get("TotalUp")->Clone();
                TH1F *theory_down=(TH1F*)pFile->Get("TotalDown")->Clone();

                theory_up->SetDirectory(0);
                theory_down->SetDirectory(0);
                pFile->Close();

                theory_up->Multiply(theory_up);
                theory_down->Multiply(theory_down);             

                local_exp_up->Add(theory_up);
                local_exp_down->Add(theory_down);

		TH1F *pHup=GetSqrtTH1(local_exp_up); //syst_up2[ff][var]);
                TH1F *pHdown=GetSqrtTH1(local_exp_down); //syst_down2[ff][var]);

		local_up.push_back(pHup);
                local_down.push_back(pHdown);

                  
                if (variable[var]->GetYield()==true){
		    std::cout <<"-- up " << pHup->Integral()  << ", down: " << pHdown->Integral() << std::endl;

	            std::vector<double> syst;
                    syst.push_back(pHup->Integral()); syst.push_back(pHdown->Integral());

                    TotalBkg[final_state[ff]]->SetSystematic(syst);
                  

                 }


          }

	  total_syst_up.push_back(local_up);
          total_syst_down.push_back(local_down);
  }



  //std::map< TString,std::map<TString, TGraphAsymmErrors*> > totalbkg;

   for (unsigned int ff=0; ff<final_state.size(); ff++){

        std::map<TString, TGraphAsymmErrors*> localmap;

        for (unsigned int var=0; var<variable.size(); var++){

           TH1F *pH=0;

           for (unsigned int process=0; process< PhysProcess.size(); process++){
              if (PhysProcess[process]->isDataSigBkg()=="isBkg"){

                if (process==0) pH=(TH1F*)PhysProcess[process]->GetHistos(ff,var)->Clone();
                else pH->Add(PhysProcess[process]->GetHistos(ff,var));
              } 
          }

           localmap[variable[var]->GetName()]=GetTGraphErrorsSyst(pH,total_syst_up[ff][var],total_syst_down[ff][var]);


         }//end var

         m_total_bkg_syst[final_state[ff]]=localmap;

   }






   return;
}

void MiniTreeAnalyzer::SetTotalBkg(){

 
   for (unsigned int ff=0; ff<final_state.size(); ff++){

         TH1F *pH=0;

        for (unsigned int var=0; var<variable.size(); var++){
           

           if (variable[var]->GetYield()==true){

                for (unsigned int process=0; process< PhysProcess.size(); process++){
                      if (PhysProcess[process]->isDataSigBkg()=="isBkg"){
                
                            if (process==0) pH=(TH1F*)PhysProcess[process]->GetHistos(ff,var)->Clone();
                            else pH->Add(PhysProcess[process]->GetHistos(ff,var));
                      }
                 }
             }// end if var
    

         }//end var


          double total_bkg_uncertainty=0;
          double total_bkg=0;

          total_bkg=pH->IntegralAndError(1,pH->GetNbinsX(),total_bkg_uncertainty);

          Yields *p_bkg_total= new Yields(pH->Integral(),total_bkg_uncertainty);

          TotalBkg[final_state[ff]]=p_bkg_total;         

   }




 return;
}
TH1F* MiniTreeAnalyzer::GetSquareTH1F(TH1F *pH){

 TH1F *pH2=0;
 pH2=(TH1F*)pH->Clone();

 pH2->Multiply(pH2,pH2);


 return pH2;
}
TGraphAsymmErrors* MiniTreeAnalyzer::GetTGraphErrorsSyst(TH1F *pStat,TH1F *up,TH1F *down){

 TGraphAsymmErrors *pTotalSystematicBkg = new TGraphAsymmErrors;


 TH1F *pHUpTotal = (TH1F*)up->Clone();
 TH1F *pHDownTotal = (TH1F*)down->Clone();

 TH1F *pH=(TH1F*)pStat->Clone();



 for (signed int bin=0; bin < pStat->GetNbinsX(); bin++){

      double stat2=pow(pStat->GetBinError(bin+1),2);

      double y=pStat->GetBinContent(bin+1);
      double x=pStat->GetBinCenter(bin+1);
      double dx=0.5*pStat->GetBinWidth(bin+1);

      double syst2_up=pow(pHUpTotal->GetBinContent(bin+1),2);
      double syst2_down=pow(pHDownTotal->GetBinContent(bin+1),2);

      double total_up=TMath::Sqrt(stat2+syst2_up);
      double total_down=TMath::Sqrt(stat2+syst2_down);


      pTotalSystematicBkg->SetPoint(bin,x,y);
      pTotalSystematicBkg->SetPointError(bin,dx,dx,total_down,total_up);

      std::cout << "-- bin: " << bin+1 << "-- x :" << x << ", y:" << y << ", total_up :"<< total_up << ", total_dn:" << total_down  << std::endl;

 }

 //std::cout << "---- Total Bkg Systematic -> "<< pHDownTotal->Integral() << " / +" << pHUpTotal->Integral() << std::endl;




  return pTotalSystematicBkg;

}
TH1F* MiniTreeAnalyzer::GetSqrtTH1(TH1F *pH){

 TH1F *pHnew=0;
 pHnew=(TH1F*)pH->Clone();
 pHnew->Reset();

 for (unsigned int i=1; i<pH->GetNbinsX(); i++){ 
        float cont=pH->GetBinContent(i);
        pHnew->SetBinContent(i,TMath::Sqrt(cont));
 }

 return pHnew;
}
std::vector<TString> MiniTreeAnalyzer::GetFullListSystematics(){

   std::vector<TString> bkg_systematics;  


   for (unsigned int process=0; process< PhysProcess.size(); process++){ 

       if (PhysProcess[process]->isDataSigBkg()=="isBkg"){ 


           std::map<TString,TString> bkg=PhysProcess[process]->GetSystematicsVariations();
           if (bkg.empty()) bkg=GetSystematicsVariations(); 


	   for (std::map<TString,TString>::iterator it=bkg.begin(); it!=bkg.end(); ++it){ 


	      if (std::find(bkg_systematics.begin(),bkg_systematics.end(),it->first) == bkg_systematics.end() ) bkg_systematics.push_back(it->first);		   
              
   	   }//end for syst

       } //end if bkg
   }//end for process


   std::cout<< "### Getting full list of systematics for bkg:" << std::endl;
   for (unsigned int i=0; i<bkg_systematics.size(); i++ ) std::cout << "--- Systematic: " << bkg_systematics[i] << std::endl;

  return bkg_systematics;
}
