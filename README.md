# controlCartography

**Public code release for:** *A Network Science Cartography of Cognitive Control System Dynamics*

**PhD Candidate Qualifying Examination**, Rutgers University Center for Molecular and Behavioral Neuroscience, 2018

**Carrisa Cocuzza** (carrisacocuzza@gmail.com), The Cole Lab (http://www.colelab.org/)

See genUsage.txt for general usage instructions and comments

## Directory Structure
Ordered by project workflow (e.g., not alphabetical)

**Directory:** partitionDemos/

- Contains demo code for implementing the Cole-Anticevic (CA) network partition to organize functional connectivity (FC) estimates (e.g., adjacency matrices) into a network community structure. Also contains demo code for assessing empirical resting-state reassignments in relation to the CA partition and applying this empirically adjusted partition to task-evoked FC data (and to be used for analyses represented in the following directories) 

**Directory:** restTaskFC_Comparisons/

- Contains MATLAB functions for computing resting-state FC and task-state FC with Pearson's Correlation Coefficient 
- Contains demo code for assessing rest-to-task changes in FC; e.g., descriptive analyses of changes to FC architecture between resting state and task-evoked state(s); as opposed to mechanistic and/or functional analyses (see below) 

**Directory:** networkDiagnosticMetrics/

- Contains demo code for diagnosing network-mechanisms of interest. Diagnostic metrics demonstrated include: global variability coefficient (GVC), between-network variability coefficient (BVC), and network flexibility (NF) 
- See networkDiagnosticMetrics/References.txt for a full reference list associated with each metric (and additional info on public repositories outside of github) 

**Directory:** cartographicMethod/

- Contains demo code for implementing our cartographic representation of cognitive control system functioning, vis-a-vis the previously demonstrated analyses and diagnostic metrics 

**Directory:** sampleResults/

- Contains sample figures and videos associated with the project herein. See sampleResultsText.txt for associated text. 
