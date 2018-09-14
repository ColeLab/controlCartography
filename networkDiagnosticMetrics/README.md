**Public code release for:** A Network Science Cartography of Cognitive Control System Dynamics  
**PhD Candidate Qualifying Examination**, Rutgers University Center for Molecular and Behavioral Neuroscience, 2018  
**Carrisa Cocuzza** (carrisacocuzza@gmail.com), The Cole Lab (http://www.colelab.org/)  

## Network diagnostic metrics
**Purpose**: use FC estimates and the network partition to diagnose variability and flexibility.

**Directory**: Contains demo code for diagnosing network-mechanisms of interest. Diagnostic metrics demonstrated include: global variability coefficient (GVC), between-network variability coefficient (BVC), and network flexibility (NF)
    - See networkDiagnosticMetrics/References.md for a full reference list associated with each metric (and additional info on public repositories outside of github)

## References & source information
1. **Global variability coefficient (GVC):** originally by Cole et al., 2013; quantifies variability in connectivity patterns across task states 
    - See References.md
    - **Further info & source code**: http://www.colelab.org/cole-etal-2013/#analysiscode
    
2. **Between-network variability coefficient (BVC):** inspired by Cole et al., 2013 and Ito et al., 2017; addresses the concern that within-network connections bias GVC results.
    - In the present study's dataset, BVC and GVC were highly comparable, suggesting that within-network connectivity does not overly bias the computation of GVC.
    - **MATLAB code:** bvcAlgorithm.m

3. **Network flexibility (NF):** originally by Bassett et al., 2011; how often a region changes network (aka module) allegiance (normalized by possible number of network assignment changes) 
    - See References.md
    - **Further info & source formulas**:  www.pnas.org/lookup/suppl/doi:10.1073/pnas.1018985108/-/DCSupplemental

4. **Network partition deviation (NPD):** novel in the present study; uses emprically-adjusted resting-state partition as a reference and assesses how often regions deviate from this pre-defined partition across task states  
    - **MATLAB code:** npdAlgorithm.m 

5. **Other network metrics**: see https://sites.google.com/site/bctnet/

## Supplementary and helper code
1. **netStruct.m:** uploads (into MATLAB workspace) useful CA partition info needed by other functions (color schemes, node indices, etc.)
    - baseDir variable needs to be edited to match your local machine's path info  

2. **restPartitionAdjuster.m:** function for adjusting CA partition by empirical resting-state preferences
    - also see the following directory: /controlCartography/restTaskFC_Comparisons
