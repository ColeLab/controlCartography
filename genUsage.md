**Public code release for:** A Network Science Cartography of Cognitive Control System Dynamics  
**PhD Candidate Qualifying Examination**, Rutgers University Center for Molecular and Behavioral Neuroscience, 2018  
**Carrisa Cocuzza** (carrisacocuzza@gmail.com), The Cole Lab (http://www.colelab.org/)  

## General Usage

### Formatting 
- All demo code files use the extension .ipynb. This corresponds to iPython Notebook (aka Jupyter Notebook). More info here: http://jupyter.org/
- All master script files are either MATLAB formatted with the .m extension, or python formatted with the .py extension. Adjacency matrices should be formatted as .mat files for MATLAB and .csv for python. 
- All demos and master scripts include code for visualization, but a directory of sample figures can be found at sampleResults/ for comparison with independent datasets. All example figures use the .png extension and are vector graphics. All videos and gifs use the .mp4 and .gif extensions, respectively.
- All demos and master scripts include statistical testing when appropriate (README.md files are included in each primary directory for further explication). 

### Project workflow & usage comments
1. **Network Partition:** input interregional functional connectivity estimates for rest and task, and obtain the empirically-adjusted Cole-Anticevic partition. 
    - Replicating this project assumes that adjacency matrices have already been estimated (to allow for different methods of fMRI pre-processing and FC estimation). The referenced work utilized standard Pearson's correlation coefficients on pre-processed BOLD activity (see manuscript for more details). While voxel- or vertex-wise data may be utilized, a parcellation scheme is advised. The referenced work utilized the 2016 Glasser parcellation: **Glasser, M. F., Coalson, T. S., Robinson, E. C., Hacker, C. D., Harwell, J., Yacoub, E., … Van Essen, D. C. (2016). A multi-modal parcellation of human cerebral cortex. Nature, 536(7615), 171–178. https://doi.org/10.1038/nature18933.** Altogether, these methods ultimately generate a 360 region x 360 region functional connectivity matrix, per participant, and per state. Our task data involved 12 task states and 50 participants (in the discovery set; another n=50 in the replication set) thus, the input of task-evoked FC matrices are of the size 360x360x50x12. Resting-state FC matrices represent 1 state per participant, and are thus of size 360x360x50. Note that all mathematical operations (averaging, standard deviation, etc.) are performed with the Fisher-Z transform for stabilizing variance across participants. 
    - We suggest that independent datasets include both resting-state and task-evoked data. However, the code herein may be adjusted for projects that only obtained task-evoked data. In the latter case we suggest using the CA partition with no adjustments. 
2. **Rest-to-Task Comparison**: compare resting state and task state FC estimates. 
    - Uses nonparametric permutation testing as well as the Mantel test. 
3. **Network Diagnostic Metrics**: use FC estimates and the network partition to diagnose variability and flexibility. 
    - See the networkDiagnosticMetrics/README.md for extensive information and source material. 
4. **Cognitive Control Cartography**: map out control control mechanisms of putative control networks, vis-a-vis the aforementioned diagnostic metrics. 
