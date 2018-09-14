**Public code release for:** A Network Science Cartography of Cognitive Control System Dynamics  
**PhD Candidate Qualifying Examination**, Rutgers University Center for Molecular and Behavioral Neuroscience, 2018  
**Carrisa Cocuzza** (carrisacocuzza@gmail.com), The Cole Lab (http://www.colelab.org/)  

## Network diagnostic metrics
**Purpose**: use FC estimates and the network partition to diagnose variability and flexibility.

**Directory**: Contains demo code for diagnosing network-mechanisms of interest. Diagnostic metrics demonstrated include: global variability coefficient (GVC), between-network variability coefficient (BVC), and network flexibility (NF)
    - See networkDiagnosticMetrics/References.txt for a full reference list associated with each metric (and additional info on public repositories outside of github)

## References & source information
1. **Global variability coefficient (GVC):** originally by Cole et al., 2013; quantifies variability in connectivity patterns across task states 
    - Cole, M. W., Reynolds, J. R., Power, J. D., Repovs, G., Anticevic, A., & Braver, T. S. (2013). Multi-task connectivity reveals flexible hubs for adaptive task control. Nature Neuroscience, 16(9), 1348–1355. https://doi.org/10.1038/nn.3470
    - **Further info & source code**: http://www.colelab.org/cole-etal-2013/#analysiscode
    
2. **Between-network variability coefficient (BVC):** inspired by Cole et al., 2013 (above reference); addresses the concern that within-network connections bias GVC results.
    - In the present study's dataset, BVC and GVC were highly comparable, suggesting that within-network connectivity does not overly bias the computation of GVC.

3. **Network flexibility (NF):** originally by Bassett et al., 2011; how often a region changes network (aka module) allegiance (normalized by possible number of network assignment changes) 
    - Bassett, D. S., Wymbs, N. F., Porter, M. A., Mucha, P. J., Carlson, J. M., & Grafton, S. T. (2011). Dynamic reconfiguration of human brain networks during learning. Proceedings of the National Academy of Sciences of the United States of America, 108(18), 7641–7646. https://doi.org/10.1073/pnas.1018985108
    - Bassett, D. S., Wymbs, N. F., Rombach, M. P., Porter, M. A., Mucha, P. J., & Grafton, S. T. (2013). Task-based core-periphery organization of human brain dynamics. PLoS Computational Biology, 9(9), e1003171. https://doi.org/10.1371/journal.pcbi.1003171
    - **Further info & source formulas**:  www.pnas.org/lookup/suppl/doi:10.1073/pnas.1018985108/-/DCSupplemental

4. **All other metrics**: see https://sites.google.com/site/bctnet/ for source code for all other network metrics utilized in the present study 
