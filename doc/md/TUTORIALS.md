BETSE Tutorial Simulations
===========

## BETSE Model Config Files

Here are several BETSE config files that will allow you to reproduce simulations
from published manuscripts, as well as become more familiar with BETSE usage.

Note that because the version of BETSE may have changed since the publication 
of the paper, it is not possible to get exactly the same results that are 
reported, however, the general concept remains the same. 

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


[Electrochemical Attractors](https://www.dropbox.com/s/m9jcon8wz8e1529/Attractors.zip?dl=0)

**Exploring Emergent Bioelectrical Patterning Mechanisms**

[Electrochemical Patterning](https://www.dropbox.com/s/0z449tn6p9c3g14/Patterns.zip?dl=0)

**Working with Electrochemical Network Simulations**

 [Electrochemical Network](https://www.dropbox.com/s/7fme976yvq2kwhw/Physiology.zip?dl=0)

## BETSE Developer's Tutorial

BETSE is a simulation software, but it is also possible to use modules and 
methods from BETSE as a toolbox in your own projects. The BETSE Developer's 
Tutorial is a Jupyter Notebook script that allows you to import BETSE modules, 
run a sim in the Notebook, and gives some pointers on how to access and plot
simulation data and data structures of the BETSE sim. 

**How to Use the Developer's Tutorial Script**

To use these simulations, first install BETSE (if you haven't done so already):

https://gitlab.com/betse/betse#simple

Download one of the available simulation packages and unzip it into a 
directory of your choice on your computer (e.g. 'BETSE_Sims'):

Link to [Developer's Tutorial](https://www.dropbox.com/s/f4ilizqnbmn2of5/developer.zip?dl=0)





