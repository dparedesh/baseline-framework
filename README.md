# Introduction
This tool was developed to perform data selection and data visualization, and it has been used to analyze the data taken by the ATLAS detector from the proton-proton collisions at the Large Hadron Collider (LHC) at CERN.

## Main settings





## 1D distributions:

Examples of plotting 1D distributions can be seen in Figure 1. By default, when plotting the signal hypothesis this one is overlaid over the background (Figure 1 (left)). However, the tool provides an option to plot the signal hypothesis stacked over the background by setting the atribute

        analyzer.doStackSignal = true  
        
in the `MiniTreeAnalyzer` (Figure 1 (middle)). When plotting real data and the expected background, it is useful to check for the compatibility of those. This can be done by computing the ratio of data over expected background by setting the option  

        analyzer.doRatioDataBkg = true

or by computing the significance of observing `n` data events given the background prediction in each bin of the histogram with

        analyzer.doSignificance = true

The resulting plot is shown at the bottom of the canvas (Figure 1 (right)). 

<p align="center">
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/1f626a01-fa09-49a6-834b-faca2d1efd1f" align="center" height="200"  ></a>
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/d4f9488a-2086-4d58-95a3-4e6ccd2daca3" align="center" height="200"  ></a>
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/0840a069-6610-491f-a4c0-954371be8123" align="center" height="190"  ></a>
</p>
<h4 align="center"><sub>Figure 1: Examples of 1D distributions created with `Tools` </sup></h4>  


## 2D distributions:

Functionality for plotting 2D distributions is also available. Examples of those distributions look like the ones below

<p align="center">
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/9b35bca8-c513-4090-952f-b31c17fa4e20" align="center" height="200"  ></a>
<a href="url"><img src="https://github.com/dparedesh/baseline-framework/assets/13987503/83a33a45-6637-4d92-bbe9-7ebff2f538a0" align="center" height="200"  ></a>
</p>
<h4 align="center"><sub>Figure 2: Examples of 2D distributions created with `Tools` </sup></h4>  


## Optimization of data selection: 

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












