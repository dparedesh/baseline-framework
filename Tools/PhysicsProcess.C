#include "PhysicsProcess.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "TMath.h"

#include "PhysicsSample.h"
#include "MiniTreeAnalyzer.h"


PhysicsProcess::PhysicsProcess():

 m_name(""),
 m_latex(""),
 m_title(""),
 m_color(0),
 m_process(""),
 m_line(0),
 m_cut(""),
 m_tree(""),
 m_weight(),
 m_positions_mass(),
 m_comb_sig(-999),
 m_is_datadriven(0),
 m_histos(),
 m_cutflow(),
 m_yield(),
 m_dict(0),
 m_Systematics(),
 m_histos_syst(),
 List_Systematics(),

 Samples()


{}
PhysicsProcess::~PhysicsProcess()
{

    for (unsigned int i=0; i< Samples.size(); i++){
        delete Samples[i];
    }
    

    for (unsigned int i=0; i< m_histos.size(); i++){
        for (unsigned int j=0; j< m_histos[i].size(); j++){
            delete m_histos[i][j];
        } // endfor  j
    } // end for i

    for (unsigned int i=0; i< m_histos_2D.size(); i++){
        for (unsigned int j=0; j< m_histos_2D[i].size(); j++){
            delete m_histos_2D[i][j];
        } // endfor  j
    } // end for i    


    for (std::map<TString,Yields*>::iterator it=m_yield.begin(); it!=m_yield.end(); ++it){
        delete it->second;
    }
    

    for (std::map<TString,std::vector< std::vector<TH1F*> > >::iterator it=m_histos_syst.begin(); it!=m_histos_syst.end(); ++it){
            
       for (unsigned int i=0; i< m_histos_syst[it->first].size(); i++){
        for (unsigned int j=0; j< m_histos_syst[it->first][i].size(); j++){
            delete m_histos_syst[it->first][i][j];
        } // endfor  j
       } // end for i*/
    }

    /*for (unsigned int i=0; i<m_cutflow.size(); i++){
        delete m_cutflow[i];
    }*/
    
    
    Samples.clear();
    m_histos.clear();
    m_yield.clear();
    m_cutflow.clear();
    m_Systematics.clear();
    m_histos_syst.clear();
    List_Systematics.clear(); 

   
} // end destructor 

void PhysicsProcess::AddSample(TString file,TString path,float scale_factor)
{


  double xs=1;


  PhysicsSample *pSample = new PhysicsSample(path+file,xs);
 
  pSample->SetScaleFactor(scale_factor);
 
  Samples.push_back(pSample);
  
  return;
 
}
void PhysicsProcess::SetStyle(TString name,TString title,TString LabelLatex,int color,int line,TString process)
{
 m_line=line;
 m_name=name;
 m_latex=LabelLatex;
 m_color=color;
 m_process=process;
 m_title=title;

 return;
}
void PhysicsProcess::SetSystematicVariations(TString list_shapes,TString list_norm){


   std::vector<TString> line_shapes;
   std::vector<TString> line_norm;

   if (list_shapes != "") line_shapes=FillVector(list_shapes);
   if (list_norm != "" ) line_norm=FillVector(list_norm);

   std::vector<TString> shapes;
   std::vector<TString> norm;

   MiniTreeAnalyzer analyzer;

   for (unsigned int i=0; i<line_shapes.size(); i++){

        std::vector<std::string> temp;
        std::string sample(line_shapes[i].Data());

        analyzer.tokenizeString(sample,'=',temp);

        shapes.push_back(temp[0]);

        List_Systematics[temp[0]]=temp[1];
   }

   for (unsigned int i=0; i<line_norm.size(); i++){

        std::vector<std::string> temp;
        std::string sample(line_norm[i].Data());

        analyzer.tokenizeString(sample,'=',temp);

        norm.push_back(temp[0]);

        List_Systematics[temp[0]]=temp[1];
   }


   for (unsigned int i=0; i<shapes.size(); i++) m_Systematics[shapes[i]]="isShape";
   for (unsigned int j=0; j<norm.size(); j++) m_Systematics[norm[j]]="isNorm";




 return;

}
void PhysicsProcess::SetSystematic(int fs,int var,TString ch,TString theory_unc){


  std::map<TString,vector< vector<TH1F*> > > LocalMap=GetMapSystematics();

   double syst2_up=0;
   double syst2_down=0;

   double yield=GetYield(ch);


   for (std::map<TString,vector< vector<TH1F*> >>::iterator it=LocalMap.begin(); it!=LocalMap.end(); ++it){ 
     
      TString syst=it->first;    
   
      double yield_syst=LocalMap[syst][fs][var]->Integral(); //GetHistosSystematics(syst,fs,var)->Integral();

      double sigma=fabs(yield_syst-yield);

      std::cout << " yield: " << yield << "-- syst: " << syst << " -- variation: " << yield_syst << std::endl;
 
      if (syst.Contains("Up") || syst.Contains("UP") || syst.Contains("up") || syst.Contains("1up" )) syst2_up=syst2_up+pow(sigma,2);
      else if (syst.Contains("Down") || syst.Contains("DOWN") || syst.Contains("down") || syst.Contains("1dn" ) ) syst2_down=syst2_down+pow(sigma,2);
      else {
              syst2_up=syst2_up+pow(sigma,2);
	      syst2_down=syst2_down+pow(sigma,2);
       }

      if (syst.Contains("JET_JER_SINGLE_NP")) syst2_down=syst2_down+pow(sigma,2);       
   }// end for syst



   TFile *pFile = new TFile(theory_unc); //input is giving in %
   TH1F *pTheory=(TH1F*)pFile->Get(ch+"_"+m_title);

   float theory_up=0,theory_down=0;

   if (pTheory){
      theory_up=pTheory->GetBinContent(1)*yield;
      theory_down=pTheory->GetBinContent(2)*yield;
   }


  std::cout<< " file :" << theory_unc << " - sample: " << m_title << std::endl;
  std::cout<<"############################## Up:" << theory_up << std::endl;

   theory_up=theory_up*theory_up;
   theory_down=theory_down*theory_down;


   syst2_up=TMath::Sqrt(syst2_up+theory_up);
   syst2_down=TMath::Sqrt(syst2_down+theory_down);

   std::vector<double> total_syst={syst2_up,syst2_down};
   

   m_yield[ch]->SetSystematic(total_syst);



 return;
}
std::vector<TString> PhysicsProcess::FillVector(TString name){

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
