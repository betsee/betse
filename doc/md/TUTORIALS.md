BETSE Tutorial Simulations
===========

## BETSE Model Config Files

The following BETSE config files that will allow you to reproduce simulations
from published manuscripts, as well as become more familiar with BETSE usage.

Note that because the version of BETSE may have changed since the publication 
of the paper, it may not be possible to get exactly the same results that are 
reported, however, the general concept remains the same in all cases. 

**How to Run Model Config Files**

To use these simulations, first install BETSE (if you haven't done so already):

https://gitlab.com/betse/betse#simple

Download one of the available simulation packages and unzip it into a 
directory of your choice on your computer (e.g. 'BETSE_Sims'):

[Electrochemical Attractors](https://www.dropbox.com/s/m9jcon8wz8e1529/Attractors.zip?dl=0)

[Electrochemical Patterning](https://www.dropbox.com/s/0z449tn6p9c3g14/Patterns.zip?dl=0)

[Electrochemical Network](https://www.dropbox.com/s/7fme976yvq2kwhw/Physiology.zip?dl=0)

Next, at the command line, switch to the directory where you have placed the 
simulation package folder with the ‘yaml’ config files and run the following 
steps:

1. Make the computational mesh:

> betse seed patterns_2018.yaml   

2. Perform the initialization:

> betse init patterns_2018.yaml 

3. Plot the initialization (optional):

> betse plot init patterns_2018.yaml 

4. Perform the simulation, which adds Nicotine from the global boundaries:

> betse sim patterns_2018.yaml 

5. Plot the simulation (optional):

> betse plot sim patterns_2018.yaml 

Repeat the above with the different config files. Note that not all of the
simulations require both an init and a sim to be run. 

## Notes on BETSE Model Config Files

**Vmem and Ion Concentrations as Attractors of Bioelectrochemical System**

This BETSE simulation demonstrates how physiological values of V<sub>mem</sub> 
and ion concentrations develop as "attractors" of the bioelectrochemical system: 
even if the system starts with completely non-physiological values (e.g.equal 
ion concentrations inside and out of the cell), it will develop 
physiologically-relevant values when it reaches steady-state.

This replicates the simulation shown in [Figure 5](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4933718/figure/F5/)
of:

> Pietak and Levin. *Exploring Instructive Physiological Signaling with the 
> Bioelectric Tissue Simulation Engine*. Front Bioeng Biotechnol. 
> 2016 Jul 6;4:55. doi: 10.3389/fbioe.2016.00055. eCollection 2016
> https://www.frontiersin.org/articles/10.3389/fbioe.2016.00055/full


(but using v1.0 of BETSE instead of v0.3)  

Download the simulation package from:

[Electrochemical Attractors](https://www.dropbox.com/s/m9jcon8wz8e1529/Attractors.zip?dl=0)

**Exploring Emergent Bioelectrical Patterning Mechanisms**

[Electrochemical Patterning](https://www.dropbox.com/s/0z449tn6p9c3g14/Patterns.zip?dl=0)

**Working with Electrochemical Network Simulations**

 [Electrochemical Network](https://www.dropbox.com/s/7fme976yvq2kwhw/Physiology.zip?dl=0)

## BETSE Developer's Tutorial

BETSE is a self-contained software package, but it is also possible to use 
modules and methods from BETSE as tools in your own projects. The BETSE Developer's 
Tutorial is a Jupyter Notebook script that allows you to import BETSE modules, 
run a sim in the Jupyter Notebook, and gives some pointers on how to access
and plot simulation data and data structures of the BETSE sim. 

**How to Use the Developer's Tutorial Script**

To use these simulations, first install BETSE (if you haven't done so already):

https://gitlab.com/betse/betse#simple

Note that using Anaconda to install BETSE will also install [Jupyter Notebooks](https://jupyter.org/). 

Download the Developers Tutorial and unzip it into a 
directory of your choice on your computer (e.g. 'BETSE_Sims'):

[Developer's Tutorial](https://www.dropbox.com/s/f4ilizqnbmn2of5/developer.zip?dl=0)

Open the Jupyter Notebook 'Dev_Demo_BETSEv1.0.ipynb' and follow instructions 
therein. Features of the BETSE simulation that the Developer's Tutorial runs
can be edited by making changes to the 'sim_config.yaml' file included with 
the Developer's Tutorial package. 





