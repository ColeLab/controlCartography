**Public code release for:** *Flexible coordinator and switcher hubs for adaptive task control.* Cocuzza et al., 2020. Journal of Neuroscience.
**Carrisa Cocuzza** (carrisacocuzza@gmail.com), The Cole Lab (http://www.colelab.org/)  

## Network Partition 
**Purpose**: input interregional functional connectivity estimates for rest and task, and obtain the empirically-adjusted Cole-Anticevic partition

**Directory**: Contains demo code for implementing the Cole-Anticevic (CA) network partition to organize functional connectivity (FC) estimates (e.g., adjacency matrices) into a network community structure. Also contains demo code for assessing empirical resting-state reassignments in relation to the CA partition and applying this empirically adjusted partition to task-evoked FC data (and to be used for analyses represented in the following directories)
  - **restPartitionAdjuster.m:** a function to apply the CA partition to FC estimates (that are parcellated via the Glasser, 2016 scheme), and adjust by empirical resting-state data preferences 
  - all helper code and variables are included in the folder 
