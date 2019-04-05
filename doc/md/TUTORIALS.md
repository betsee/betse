BETSE Tutorial Simulations
===========

## BETSE Model Config Files

The following BETSE config files allow you to reproduce simulations
from published manuscripts, as well as become more familiar with BETSE usage.

Note that because the version of BETSE may have changed since the publication 
of the paper, it may not be possible to get exactly the same results that are 
reported in the paper, however, the general concept remains the same in all cases. 

**How to Run Model Config Files**

To use these simulations, first install BETSE (if you haven't done so already):

https://gitlab.com/betse/betse#simple

Download one of the available simulation packages and unzip it into a 
directory of your choice on your computer:

[Electrochemical Attractors](https://www.dropbox.com/s/m9jcon8wz8e1529/Attractors.zip?dl=0)

[Electrochemical Patterning](https://www.dropbox.com/s/0z449tn6p9c3g14/Patterns.zip?dl=0)

[Electrochemical Network](https://www.dropbox.com/s/7fme976yvq2kwhw/Physiology.zip?dl=0)

Next, at the command line, switch to the directory where you have placed the 
simulation package folder with the ‘yaml’ config files and run the following steps:

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
simulations require both an 'init' and a 'sim' to be run (as discussed below). 

## Notes on BETSE Model Config Files

**Vmem and Ion Concentrations as Attractors of Bioelectrochemical System**

This BETSE simulation demonstrates how physiological values of V<sub>mem</sub> 
and ion concentrations develop as "attractors" of the bioelectrochemical system. 
Even if the electrochemical cell system starts with completely 
non-physiological values (e.g.equal ion concentrations inside and out of the 
cell), it will develop physiologically-relevant values when it reaches 
steady-state.

This replicates the simulation shown in [Figure 5](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4933718/figure/F5/)
of:

> Pietak and Levin. *Exploring Instructive Physiological Signaling with the 
> Bioelectric Tissue Simulation Engine*. Front Bioeng Biotechnol. 
> 2016 Jul 6;4:55. doi: 10.3389/fbioe.2016.00055. eCollection 2016
> https://www.frontiersin.org/articles/10.3389/fbioe.2016.00055/full


*(originally using BETSE version 0.3)*  

Download the simulation package from:

[Electrochemical Attractors](https://www.dropbox.com/s/m9jcon8wz8e1529/Attractors.zip?dl=0)

**Exploring Emergent Bioelectrical Patterning Mechanisms**

Interesting patterns linking bioelectrical and molecular states can self-assemble
when a "gating electrodiffusion" mechanism is active. Here, a charged molecule 
gates an ion channel to alter V<sub>mem</sub> and also moves by electrodiffusion between
gap junction coupled cells of the cluster in the V<sub>mem</sub> gradients it
generates. This produces a range of categorically different Vmem patterns. 
Some of the results of this model are presented in [Figure 4](https://www.sciencedirect.com/science/article/pii/S0079610718300415?via%3Dihub#fig4)
of the paper:

> Pietak and Levin. Bioelectrical control of positional information in development 
> and regeneration: A review of conceptual and computational advances. Progress in 
> Biophysics and Molecular Biology. 2018. 137:52-68.
> https://www.sciencedirect.com/science/article/pii/S0079610718300415?via%3Dihub


*(originally using BETSE version 0.7)*

To obtain different patterns, there are only a small number of variables that 
alter pattern, importantly, many of these are in the 'worm_3.yaml' file of the 
'extra_configs' directory, included with this simulation package. In the 
simulation package config files, parameters that are worth altering to change 
patterns are marked by an '@' symbol, with a useful range given. 

Note: There is no real need to run the "sim" phase for these -- the system is 
set up so that it should form a pattern using only the 'init' phase command. 

Download the simulation package from:

[Electrochemical Patterning](https://www.dropbox.com/s/0z449tn6p9c3g14/Patterns.zip?dl=0)

**Working with Electrochemical Network Simulations**

This simulation features a small circular cluster of cells simulating numerous 
interventions applied to the cells for 5 minute intervals to show the effect of
these activities on V<sub>mem</sub>, pH, and metabolic status.  

This config file presents a number of examples for how to do things like change 
the K+ level of the environment, dynamically change concentrations of custom 
substances at the environment, and define and work with your own custom 
ion pumps/transporters and channels. 

During the 'sim' phase, the following interventions are applied:
* 60 -120 seconds, drug 'X' is added to the environment and blocks a K+ leak channel, cells depolarize 
* 240 -300 seconds, drug 'Y', a Na+ channel agonist, is added to the environment, cells depolarize
* 420 - 480 seconds, the Na/K-ATPase pump is blocked by addition of ouabain, cells depolarize slightly
* 600 - 660 seconds, cAMP is added to the environment, which activates an H,K-ATPase pump, cell pH increases
* 780-840 seconds, a potassium salt, K+/M- is added to the environment, cells depolarize

all agents are added to the environment by transiently changing the environmental boundary
condition for their concentration.

Results from this model are presented in [Figure 3](https://www.sciencedirect.com/science/article/pii/S0079610718300415?via%3Dihub#fig3)
of the paper:

> Pietak and Levin. Bioelectrical control of positional information in development 
> and regeneration: A review of conceptual and computational advances. Progress in 
> Biophysics and Molecular Biology. 2018. 137:52-68.
> https://www.sciencedirect.com/science/article/pii/S0079610718300415?via%3Dihub

*(originally using BETSE version 0.7)*

Download the simulation package from:

[Electrochemical Network](https://www.dropbox.com/s/7fme976yvq2kwhw/Physiology.zip?dl=0)
 

