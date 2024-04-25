## Introduction
This tool was developed to perform data selection and data visualization, and it has been used to analyze the data taken by the ATLAS detector from the proton-proton collisions at the Large Hadron Collider (LHC) at CERN. The tool has been developed on top of the [ROOT Data Analysis Framework](https://root.cern/), a high-performance software written mainly in C++. 

## Description

The tool consists of several classes. The ones at the forefront of the user are:

- `Channel`:  It describes the selection to be applied to the samples.
- `VariableDistr`: The variables used to create the plots.
- `MiniTreeAnalyzer`:  This is the main class, and it is the one used to execute the tool. It uses as input the `channels` and `variables` to plot, as well as the samples to be used. All main settings are done in this class.

After executing the tool you get access to additional objects created during execution: 

- `PhysicsProcess`: It corresponds to one or more samples that are summed up and treated equally when adding the process to the `MiniTreeAnalyzer` (see later). 
- `Yields`: It holds the number of events for every physics process, together with its statistical and systematic uncertainty. 

## Setting the tool:

### Creating a `Channel`:

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


### Creating a `VariableDistr`:

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

  - `name`: File containing the list of samples (file format must be in .root). Examples of the lists to be provided can be found [here](https://github.com/dparedesh/baseline-framework/tree/master/InputList)
  - `title`: Name of the object shown at the terminal when running.
  - `latex`: Custom name shown in the plot legend.
  - `color`: Integer referencing the colors in [TColor](https://root.cern.ch/doc/master/classTColor.html)
  - `process`: It can have three different values: `"isSig"`, `"isBkg"`, and `"isData"`. Depending on the option provided the physics process will be plotted in a different style, i.e.
 
    - `isSig`: The physics process will be shown as a line overlaid to the background distribution.
    - `isBkg`: The physics process will be stacked for all background processes. Histograms will be color-filled.
    - `isData`: Data will be shown as black dot points with their bar uncertainty.  
  - `isDD`: Samples of simulated events are usually weighted. If they are estimated using data-driven techniques (`isDD`), then the weight must not be applied.
  - `cut`: If provided, the selection applied to this physics process will be given by this variable and not by the one defining the channel.
  - `weight`: If provided, the weight applied to this physics process will be given by this variable and not by the one defined in the analyzer (see later).
  - `tree`:  If provided, the TTree where this physics process is stored will be given by this variable and not by the one defined in the analyzer.
 
- Additional settings: Depending on what you want to do, you will need to set extra options in the `MiniTreeAnalyzer`. The most usual are:

  - Weight to be applied to all samples of simulated events:

    ```cpp
    
    AddWeight(TString weight)
    ```
  - By default, the tool is used to print the number of events per variable and per channel. If you wish to show the distributions you must tell it to the analyzer via the attribute

    ```cpp
    
    doPlotDistributions = true
    ```
  - Plots can be saved if requested. This is done with the following function:
 
    ```cpp
    
    SavePlots(bool save,TString folder)
    ```
    where `folder` is the directory where the plots will be saved.

  - A log version of the plots can be performed if requested. It must be done via the attribute

    ```cpp
   
    doLogPlots = true

    ```  
  - Histograms in ROOT format can be saved if requested via the function:

    ```cpp
 
    SaveHistos(bool save,TString folder="Histos")
    ```     
    where `folder` is the name of the directory where the .root files will be saved.

  - ATLAS Collaboration requires labeling the plots according to the approval level in the collaboration. It also requires to show the integrated luminosity used for the data sample. This is done with the following functions:
 
    ```cpp

    SetATLASLabel("Internal");
    SetLuminosity(139,"139 fb^{-1}");    
    ```
  - For debugging and print more information, it is useful to set the following function:

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

### 1D distributions:


#### Default output:

Examples of the default 1D distributions can be seen in Figure 1.   By default, when plotting the signal hypothesis this one is overlaid over the background. 

<p align="center">
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/1f626a01-fa09-49a6-834b-faca2d1efd1f" align="center" height="210"  ></a>
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/273f3723-bb4d-4a4c-a1c5-02fbdf8fbd39" align="center" height="210"  ></a>
</p>
<h4 align="center"><sub>Figure 1: Examples of the default 1D distributions created with `Tools`. The signal is represented by the blue dotted line.  </sup></h4>  

 <br/><br/> 


#### Stacked signal:
 
 The tool provides an option to plot the signal hypothesis stacked over the background by setting the attribute

 ```cpp        
 analyzer.doStackSignal = true  
 ```

in the `MiniTreeAnalyzer`. An example showing the result can be seen in Figure 2. 

<p align="center">
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/a37d261b-f24d-45e4-a36e-f1ee8a46279d" align="center" height="190"  ></a>
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/d4f9488a-2086-4d58-95a3-4e6ccd2daca3" align="center" height="200"  ></a>
</p>
<h4 align="center"><sub>Figure 2: Examples of 1D distributions created with `Tools` using stacked signal hypothesis. The signal is represented by shaded histograms stacked on top of the total background. </sup></h4>  


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
<h4 align="center"><sub>Figure 3: Examples of 1D distributions created with `Tools` showing the ratio of the  data and the total background. </sup></h4>  


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
<h4 align="center"><sub>Figure 4: Examples of 1D distributions created with `Tools` showing the significance of observing `n` data events given the background prediction in each bin of the histogram. </sup></h4>  


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
<h4 align="center"><sub>Figure 5: Examples of 1D distributions showing the shape comparison of the different physics processes. Distributions are normalized to unit.  </sup></h4>  


<br/><br/>

### 2D distributions:

Functionality for plotting 2D distributions is also available. It can be obtained the passing the variables to be plotted to the following function of the analyzer:

```cpp
 AddVariable2D(VariableDistr* var_x,VariableDistr* var_y)
 ```

Examples of 2D distributions obtained with the tool are shown in Figure 6.

<p align="center">
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/9b35bca8-c513-4090-952f-b31c17fa4e20" align="center" height="220"  ></a>
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/83a33a45-6637-4d92-bbe9-7ebff2f538a0" align="center" height="220"  ></a>
</p>
<h4 align="center"><sub>Figure 6: Examples of 2D distributions created with `Tools` </sup></h4>  

<br/><br/>

### Tables

If requested, the tool will provide the yield tables in Latex format. Those yiels will include not only the nominal counting but also the statistical and systematic uncertainties for every physics process and the total background. This can be done by calling the following function when setting the analyzer:

```cpp
 SaveYieldsTables(bool save,TString folder="Tables")
```
where `folder` indicates the directory where the tables will be saved.  The table produced looks like the one shown in Figure 7.

<p align="center">
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/f996f866-c13d-4270-a743-cac5c9692715" align="center" height="220"  ></a>
</p>
<h4 align="center"><sub>Figure 7: Table generated with `Tools` containing the nominal event counting and its statistical and systematic uncertainty. Values are computed for every  process and for every decay channel (data selection). The total expected background is also computed. The number of observed events in data is also shown. </sup></h4>  





<br/><br/>


### Optimization of data selection: 

A script can be also created inputing `strings` in a template macro as the example given here. The significance  computed as , can be used to find the best selection that maximize the separation between the signal produced by a new particle and the one below.  The results obtained with this selection can be used to perform a statitical interpretation of the data, by computing the 95% CL upper limits on the production rate of the hypothetical particle. Examples of those limits can be found below.  

<p align="center">
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/76ee9a09-dd95-4526-ac4f-22bf5dd7d276" align="center" height="200"  ></a>
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/596802b8-39f8-4197-8ef8-7c40f8103ce6" align="center" height="200"  ></a>
</p>
<h4 align="center"><sub>Figure 3: Upper limits computed at 95% CL on the production rate of the hypothetical particle. The data selected with `Tools` have been used to perform the statistical interpretation</sup></h4>  




## Official results obtained with this tool..

The tool has been used to produce the official results of the following paper publications:

- [JHEP 02 (2024) 107](https://inspirehep.net/literature/2673888)
- [JHEP 11 (2023) 150 ](https://link.springer.com/article/10.1007/JHEP11(2023)150)
- [Phys. Rev. D 100 (2019) 012006](https://inspirehep.net/literature/1711261)
- [JHEP 06 (2018) 166](https://link.springer.com/article/10.1007/JHEP06(2018)166)
- [Phys. Rev. D 92 (2015) 092001](https://inspirehep.net/literature/1373520)












