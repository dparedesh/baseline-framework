

## Introduction
This tool was developed to perform data selection and data visualization, and it has been used to analyze the data taken by the *ATLAS* detector from the proton-proton collisions at the *Large Hadron Collider (LHC)* at *CERN*. The tool has been developed on top of the [ROOT Data Analysis Framework](https://root.cern/), a high-performance software written mainly in C++. 

## Description

The tool consists of several classes. All of them can be found in the [Tools](https://github.com/dparedesh/baseline-framework/tree/master/Tools) directory of this repository. The ones at the forefront of the user are:

- `Channel`:  It describes the selection to be applied to the samples.
- `VariableDistr`: The variables used to create the plots.
- `MiniTreeAnalyzer`:  This is the main class, and it is the one used to execute the tool. It uses as input the `channels` and `variables` to plot, as well as the samples to be used. All main settings are done in this class.

After executing the tool you get access to additional objects created during execution: 

- `PhysicsProcess`: It corresponds to one or more samples that are summed up and treated equally when adding the process to the `MiniTreeAnalyzer` (see later). 
- `Yields`: It holds the number of events for every physics process, together with its statistical and systematic uncertainty. 

## Setting the tool

### Creating a `Channel`

A channel is created using the following syntax:

```cpp
Channel(TString name, TString latex,const TString cut)
```

where

- `name`: Internal identifier for the channel.
- `latex`: Custom name that will be shown as the label in the plot. Usually provided in Latex format. 
- `cut`:  Selection to be applied to the data. 

Example: 

```cpp
TString Rpc2L0b = "(nSigLep>=2 && nBJets20==0 && nJets40>=6 && met>200000 && meff>1000000 && met/meff>0.2)";
        
Channel *ch_Rpc2L0b = new Channel("Rpc2L0b","Rpc2L0b",Rpc2L0b);
```


### Creating a `VariableDistr`

A variable is created with the following syntax:

```cpp
VariableDistr(const TString &variable,const TString &title,const TString x_label, const TString y_label,
              const Int_t nbBins,const float xMin,const float xMax,
              bool yield=false,bool logX=false,bool logY=false);
```

where

- `variable`: Internal identifier.
- `title`: Custom name that will be shown at the terminal when executing the tool.
- `x_label`, `y_label`: Custom labels that will be shown at the plot for the x-axis and y-axis, respectively. To be provided in Latex format.
- `nbBins`: Number of bins of the distribution.
- `xMin`, `xMax`: The minimum and maximum values for which the data will be plotted.
- `yield`: Use this variable to compute the number of events contained in the distribution. 
- `logX`, `logY`: Uses log scale for the x-axis or y-axis, respectively.  


Example: 

```cpp

VariableDistr *var_meff = new VariableDistr("meff","meff","m_{eff} [GeV]","Events/",11,400,1500,true);
```

### Setting the `MiniTreeAnalyzer`

Once you have created the channels and the variables to be plotted, you can set the MiniTreeAnalyzer. The object is created like this:

```cpp

MiniTreeAnalyzer analyzer;
```

The `MiniTreeAnalyzer` needs a minimum input to run:

- `datDir`: Attribute of the analyzer. Directory where the samples of real data are stored.
- `bkgDir`: Attribute of the analyzer. Directory where the samples of simulated background events are stored.
- `sigDir`: Attribute of the analyzer. Directory where the samples of simulated signal events are stored. This argument is actually optional, it must be provided only if signal samples are plotted.
- The name of the TTree object containing the events in the ROOT file. This is set with the function:
  
  ```cpp

  SetTreeName(TString tree) 
  ```
 
- The `Channel` previously created. You can add sequentially as many channels as you want. This is done with the function:

  ```cpp

  AddChannel(Channel *p_ch)
  ```
- The `VariableDistr` previously created. You can add sequentially as many variables as you want with the function:

  ```cpp
  AddVariable(VariableDistr* var)
  ```
- The physics processes to be studied. You can add as many processes as you want using the following syntax:
  
  ```cpp
  
  AddProcess(TString name,TString title,TString latex,int color,int line,TString process,
             bool isDD=false,TString cut="",TString weight="",TString tree="")
  ```
  where

  - `name`: File containing the list of samples (the samples must be in the format .root). Examples of the lists to be provided can be found [here](https://github.com/dparedesh/baseline-framework/tree/master/InputList).
  - `title`: Name of the object shown at the terminal when running.
  - `latex`: Custom name shown in the plot legend.
  - `color`: Integer referencing the colors in [TColor](https://root.cern.ch/doc/master/classTColor.html)
  - `process`: It can have three different values: `"isSig"`, `"isBkg"`, or `"isData"`. Depending on the option provided the physics process will be plotted in a different style, i.e.
 
    - `isSig`: The physics process will be shown as a line overlaid to the background distribution.
    - `isBkg`: The physics process will be stacked for all background processes. Histograms will be color-filled.
    - `isData`: Data will be shown as black dot points with their bar uncertainty.  
  - `isDD`: Samples of simulated events are usually weighted. If they are estimated using data-driven techniques (`isDD`), then the weight must not be applied.
  - `cut`: If provided, the selection applied to this physics process will be given by this variable and not by the one defining the channel.
  - `weight`: If provided, the weight applied to this physics process will be given by this variable and not by the one defined in the analyzer (see later).
  - `tree`:  If provided, the TTree where this physics process is stored will be given by this variable and not by the one defined in the analyzer.
 
- Additional settings: Depending on what you want to do, you will need to set extra options in the `MiniTreeAnalyzer`. The most usual are:

  - Apply a weight to the samples: In general, all samples generated via Montecarlo simulation must be weighted accordingly. This can be done via the function:

    ```cpp
    
    AddWeight(TString weight)
    ```
  - Plot distributions: By default, the tool prints the number of events per variable and channel. If you wish to *show* the distributions you must tell it to the analyzer via the attribute

    ```cpp
    
    doPlotDistributions = true
    ```
  - Save plots: Plots can be saved if requested. This is done with the following function:
 
    ```cpp
    
    SavePlots(bool save,TString folder)
    ```
    where `folder` is the directory where the plots will be saved. Plots will be saved in different formats. 

  - Log plots: A log version of the plots can be done if requested. It must be done via the attribute

    ```cpp
   
    doLogPlots = true

    ```  
  - Save histograms: Histograms in ROOT format can be saved if requested via the function:

    ```cpp
 
    SaveHistos(bool save,TString folder="Histos")
    ```     
    where `folder` is the name of the directory where the .root files will be saved.

  - Show *ATLAS label and luminosity*: *ATLAS Collaboration* requires labeling the plots according to the approval level in the collaboration. It also requires to show the integrated luminosity used for the data sample. This is done with the following functions:
 
    ```cpp

    SetATLASLabel("Internal");
    SetLuminosity(139,"139 fb^{-1}");    
    ```
  - Print debug info: For debugging and printing more information, it is useful to set the following function:

    ```cpp     
    SetDebugLevel(int debug=0) 
    ```

- Order the tool to start the execution:

 ```cpp 
 Execute()
 ```

Example: 

```cpp
     //Defining main settings
     MiniTreeAnalyzer analyzer;
     analyzer.datDir="/afs/cern.ch/user/d/dparedes/WorkCERN/Analysis_SUSY_SS3L_2018/HistFitterComplexFinal/SS3Lep/SS3L_HF/FullRun2/prepare/"; 
     analyzer.bkgDir="/eos/user/d/dparedes/SUSYComplex/";
     analyzer.sigDir="/eos/user/d/dparedes/SUSYComplex/";
     analyzer.SetTreeName("evtel");
     analyzer.SetATLASLabel("Internal"); 
     analyzer.SetLuminosity(139,"139 fb^{-1}");   
     analyzer.SetDebugLevel(1); 

     //Adding decay channels     
     analyzer.AddChannel(ch_Rpc2L0b);
     analyzer.AddChannel(ch_Rpc2L1b);
     analyzer.AddChannel(ch_Rpc2L2b);
     analyzer.AddChannel(ch_Rpc3LSS1b);
     analyzer.AddChannel(ch_Rpv2L0b);

     //adding variables to the analyzer
     analyzer.AddVariable(var_nlep);
     analyzer.AddVariable(var_njets);
     analyzer.AddVariable(var_meff);

     //Adding background processes
     analyzer.AddProcess("InputList/Bkg.txt","OtherMultiboson","OtherMultiboson",kGreen+2,1,"isBkg",false,"","","OtherMultiboson_nom",3);
     analyzer.AddProcess("InputList/Bkg.txt","ttH","ttH",kOrange,1,"isBkg",false,"","","ttH_nom",3); 
     analyzer.AddProcess("InputList/Bkg.txt","RareTop_NottH","RareTop (no ttH)",kCyan+1,1,"isBkg",false,"","","RareTop_NottH_nom",3);
             

     //Adding real data
     analyzer.AddProcess("InputList/Data.txt","Data","Data",1,1,"isData",false,"","","data_nom",3);

     //Adding signal generated by a new particle
     analyzer.AddProcess("InputList/signal.txt","Signal","Signal",kRed,2,"isSig",false,"","","",3);


     //Montecarlo simulation samples will be weighted by this weight
     analyzer.AddWeight("totweight*lumiScaling");   

     //More settings... 
     analyzer.SetPrecision(2);
     analyzer.SavePlots(true,"Plots");
     analyzer.doLogPlots=false;
     analyzer.SaveHistos(true,"Histos");
    
     // Finally: execute the analyzer to get the final plots and tables
     analyzer.Execute();
 
```
        

## Output of the tool and more custom options

### 1D distributions


#### Default output

Examples of the default 1D distributions can be seen in Figure 1.   By default, when plotting the signal hypothesis this one is overlaid over the background. 

<p align="center">
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/1f626a01-fa09-49a6-834b-faca2d1efd1f" align="center" height="210"  ></a>
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/273f3723-bb4d-4a4c-a1c5-02fbdf8fbd39" align="center" height="210"  ></a>
</p>
<h4 align="center"><sub>Figure 1: Examples of the default 1D distributions created with `Tools`. The signal is represented by the blue dotted line. Source [7].  </sup></h4>  

 <br/><br/> 

#### Stacked signal
 
 The tool provides an option to plot the signal hypothesis stacked over the background by setting the attribute

 ```cpp        
 analyzer.doStackSignal = true  
 ```

in the `MiniTreeAnalyzer`. An example showing the result can be seen in Figure 2. 

<p align="center">
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/a37d261b-f24d-45e4-a36e-f1ee8a46279d" align="center" height="190"  ></a>
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/d4f9488a-2086-4d58-95a3-4e6ccd2daca3" align="center" height="200"  ></a>
</p>
<h4 align="center"><sub>Figure 2: Examples of 1D distributions created with `Tools` using stacked signal hypothesis. The signal is represented by shaded histograms stacked on top of the total background. Source [6]. </sup></h4>  


<br/><br/> 

#### Showing the ratio of data and background prediction

When plotting real data and the expected background, it is useful to check for their compatibility. This can be done by computing the ratio of data over the expected background by setting the option  

 ```cpp
analyzer.doRatioDataBkg = true
 ```
The resulting plot is shown at the bottom of the canvas in Figure 3.  

<p align="center">
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/3aeff934-ba9a-4d0a-912f-631d23396d90" align="center" height="230"  ></a>
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/144ecb0e-fa0e-4f27-8a68-de618f01d632" align="center" height="230"  ></a>
</p>
<h4 align="center"><sub>Figure 3: Examples of 1D distributions created with `Tools` showing the ratio of the  data and the total background. Source [5]. </sup></h4>  


<br/><br/> 

#### Showing the significance

A more elegant way to check for the compatibility between the data and the total background can be obtained by computing the **significance** of observing `n` data events given the background prediction in each bin of the histogram. This can be done by setting the following option in the analyzer:

 ```cpp
 analyzer.doSignificance = true
 ```

The resulting plot is shown in Figure 4. 
<p align="center">
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/0840a069-6610-491f-a4c0-954371be8123" align="center" height="200"  ></a>
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/3506e95d-1d61-47f7-9c27-ff7ba7e92b05" align="center" height="200"  ></a>
</p>
<h4 align="center"><sub>Figure 4: Examples of 1D distributions created with `Tools` showing the significance of observing `n` data events given the background prediction in each bin of the histogram. Source [7]. </sup></h4>  


<br/><br/>


#### Comparing shapes

One of the most common things when comparing different physics processes is comparing their shape, while keeping the histograms normalized to unit. This can be done by setting the option 

 ```cpp
 analyzer.BuildShapesNormalized()  //for histos normalized to unit
 analyzer.BuildShapes()            //if no normalization needs to be applied
 ```

The final plots will look like the ones shown in Figure 5. 

<br/><br/>

<p align="center">
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/36015733-235e-4a1f-b921-04f2ac1cbb11" align="center" height="210"  ></a>
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/4840fcf6-9bf3-4be1-8fc2-1fb940699959" align="center" height="210"  ></a>
</p>
<h4 align="center"><sub>Figure 5: Examples of 1D distributions showing the shape comparison of the different physics processes. Distributions are normalized to unit. Source [2]. </sup></h4>  


<br/><br/>

### 2D distributions

Functionality for plotting 2D distributions is also available. It can be obtained the passing the variables to be plotted to the following function of the analyzer:

```cpp
 AddVariable2D(VariableDistr* var_x,VariableDistr* var_y)
 ```

Examples of 2D distributions obtained with the tool are shown in Figure 6.

<p align="center">
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/9b35bca8-c513-4090-952f-b31c17fa4e20" align="center" height="220"  ></a>
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/83a33a45-6637-4d92-bbe9-7ebff2f538a0" align="center" height="220"  ></a>
</p>
<h4 align="center"><sub>Figure 6: Examples of 2D distributions created with `Tools`. Source [7]. </sup></h4>  

<br/><br/>

### Tables

If requested, the tool will provide the yield tables in Latex format. Those yiels will include not only the nominal counting but also the statistical and systematic uncertainties for every physics process and the total background. This can be done by calling the following function when setting the analyzer:

```cpp
 SaveYieldsTables(bool save,TString folder="Tables")
```
where `folder` indicates the directory where the tables will be saved.  The table produced looks like the one shown in Figure 7.

<p align="center">
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/f996f866-c13d-4270-a743-cac5c9692715" align="center" height="180"  ></a>
</p>
<h4 align="center"><sub>Figure 7: Table generated with `Tools` containing the nominal event counting and its statistical and systematic uncertainty. Values are computed for every  process and for every decay channel (data selection). The total expected background is also computed. The number of observed events in data is also shown. </sup></h4>  





<br/><br/>


### Optimization of data selection using the signal significance 


One of the standard methods to find the optimal selection to be applied to the samples is maximizing the signal significance. To compute the significance of every signal hypothesis for every channel you must call the following function in the analyzer before the execution:

```cpp
SetSignificance(double rel_syst=0.3,bool root_sig=true,bool stat_forum=false)
```

where 

- `rel_syst:` Relative value of the systematic uncertainty to be applied to the total background. 
- `root_sig:` Significance as computed by ROOT with the function `RooStats::NumberCountingUtils::BinomialExpZ()`. This is the default value as it allows a fast execution of the program.
- `stat_forum:` This is the  recommended way to compute the significance. However, it is computationally more expensive. The significance is computed as recommended by the *ATLAS Statistical Forum*  in [Formulae for Estimating Significance](https://cds.cern.ch/record/2736148?ln=es):


<p align="center">
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/25e02a50-7189-4a98-8c49-b923200ff8bc" align="center" height="90"  ></a>
</p>

   where $n$ is the number of signal + background events, $b$ is the total background, and $\sigma$ is the total uncertainty on the background. 
  
 <br/><br/>
 After the execution step, you will have access to the `PhysicsProcess` objects created during the execution. You can get the significance for every signal hypothesis and for every channel by using the following code: 

```cpp
vector<PhysicsProcess*> Processes=analyzer.GetPhysicsProcesses();
vector<Channel*> Channels=analyzer.GetChannels();
std::map<TString,Yields*> Bkg=analyzer.TotalBkg;
std::map< TString,std::map<TString,double> > v_Significances=analyzer.RootSignificance;


for (unsigned int j=0; j<Channels.size(); j++){
        std::cout << "--- In channel : " << Channels[j]->GetName() << std::endl;
        PrintOutput(Channels[j],Processes,Bkg,v_Significances);
}
```

where the function `PrintOutput(<args>)` can be defined as:

```cpp

void PrintOutput(Channel *Channels,vector<PhysicsProcess*> Processes,
                 std::map<TString,Yields*> Bkg,
                 std::map<TString,std::map<TString,double> > v_Significances){
     
     for (unsigned int i=0; i<Processes.size(); i++){

          std::cout << "-- Process :" << Processes[i]->GetTitle() <<
                       ",  yield: "   << Processes[i]->GetYield(Channels->GetName()) <<
                       "+/-"          << Processes[i]->GetStatistical(Channels->GetName())  << std::endl;

          if (Processes[i]->isDataSigBkg()=="isSig"){
             std::cout << "..... Significance: " << v_Significances[Processes[i]->GetName()][Channels->GetName()] << std::endl;
          }
     } 

     return;
}

```


At this point, it is important to remember that what makes a channel different from another is the *selection applied to the samples*. Thus, different selections will output different significances. The goal is to find the selection that maximizes the significance.  

You can create different channels by scanning the values of the variables of interest and building a channel for every different selection. A script that loops over different channels and inputs that channel to a template macro that computes the significance can be used for this scan. The selection that maximizes the significance is the one that will be used to provide the final results. The statistical interpretation is then performed for the selected data by computing the 95% CL upper limits on the production rate of the hypothetical particle. Examples of those upper limits can be seen in Figure 8.  

<p align="center">
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/76ee9a09-dd95-4526-ac4f-22bf5dd7d276" align="center" height="200"  ></a>
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/596802b8-39f8-4197-8ef8-7c40f8103ce6" align="center" height="200"  ></a>
</p>

<p align="center">
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/1ca9275d-2c23-4309-b0a3-89a32b0e874e" align="center" height="220"  ></a>
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/837b2115-7362-498c-8cd3-0b06f0ac7dc2" align="center" height="220"  ></a>      
</p>
<h4 align="center"><sub>Figure 8: Upper limits computed at 95% CL on the production rate of the hypothetical particle. The data selection obtained with `Tools` has been used to perform the statistical interpretation. Source [1,6]. </sup></h4>  


<br/><br/>

## Official results obtained with this tool...

The tool has been used to produce the official results of the following paper publications:

- [[1] JHEP 02 (2024) 107](https://inspirehep.net/literature/2673888): Search for pair production of squarks or gluinos decaying via sleptons or weak bosons in final states with two same-sign or three leptons with the ATLAS detector.
- [[2] JHEP 07 (2023) 203](https://link.springer.com/article/10.1007/JHEP07(2023)203): Search for $t\bar tH/A \rightarrow t\bar tt\bar t$ production in the multilepton final state in proton-proton collisions at $\sqrt s =13$ TeV with the ATLAS detector.
- [[3] JHEP 11 (2023) 150 ](https://link.springer.com/article/10.1007/JHEP11(2023)150): Search for direct production of winos and higgsinos in events with two same-charge leptons or three leptons in pp collision data at $\sqrt s =13$ TeV with the ATLAS detector.
- [[4] JHEP 06 (2020) 046](https://link.springer.com/article/10.1007/JHEP06(2020)046): Search for squarks and gluinos in final states with same-sign leptons and jets using 139 $fb^{−1}$ of data collected with the ATLAS detector.
- [[5] Phys. Rev. D 100 (2019) 012006](https://inspirehep.net/literature/1711261): Search for chargino and neutralino production in final states with a Higgs boson and missing transverse momentum at $\sqrt s =13$ TeV with the ATLAS detector.
- [[6] JHEP 06 (2018) 166](https://link.springer.com/article/10.1007/JHEP06(2018)166): Search for Higgs boson decays to beyond-the-Standard-Model light bosons in four-lepton events with the ATLAS detector $\sqrt s =13$ TeV.
- [[7] Phys. Rev. D 92 (2015) 092001](https://inspirehep.net/literature/1373520): Search for new light gauge bosons in Higgs boson decays to four-lepton final states in $pp$ collisions at $\sqrt s =13$ TeV with the ATLAS detector at the LHC.












