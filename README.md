# NMQA-QSLAM Research Codebase

Noise Mapping for Quantum Architectures (NMQA) using Quantum Simultaneous Localisation and Mapping (QSLAM).

References:

[1] Gupta, R et. al. Adaptive scheduling of noise characterization in quantum computers (2019). Accessed https://arxiv.org/abs/1904.07225 
[2] Gupta, R & Biercuk, M. J.  Convergence analysis for NMQA (forthcoming 2019)
[3] Gupta, R. & Biercuk, M. J. Spatial sampling analysis for NMQA using Padua points (forthcoming 2019)

Contains the following packages:

    clfanalysis : Python / Artemis scripts to implement a classifier to assign single qubit binary measurements to images of trapped ions. Reference [1].
    expt_data : Process trapped ion images into measurement database. Codebase stored as a series of Jupyter Notebooks. Reference [1].
    expt_qslam : Python / Artemis scripts to run NMQA/QSLAM and compared with Naive Approach using experimental measurement databased. Reference [1].
    intrarun : Python / Artemis scripts to implement convergence analysis for NMQA-QSLAM. Reference [2].
    paduaq : Python package for Lagrange polynomial interpolation on 2D square at Padua points. Reference [3].
    paperdata : Python / Artemis scripts to compare NMQA with Naive approaches. Reference [1].
    qslam : Python package to implement NMQA algorithm. Reference [1], [2], [3].
    sampling_analysis: Python / Artemis scripts to compare NMQA with Padua interpolation techniques. Reference [3].
    statedetectn : Python / Artemis scripts for ion state detection. Reference [1].
