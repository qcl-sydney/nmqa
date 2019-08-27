# NMQA-QSLAM 

Research Codebase Noise Mapping for Quantum Architectures (NMQA). Formerly known as: Quantum Simultaneous Localisation and Mapping (QSLAM).

References:

[1] Gupta, R et. al. Adaptive scheduling of noise characterization in quantum computers (2019). Accessed https://arxiv.org/abs/1904.07225 

[2] Gupta, R & Biercuk, M. J.  Convergence analysis for NMQA (unpublished 2019)

[3] Gupta, R. & Biercuk, M. J. Spatial sampling analysis for NMQA using Padua points (unpublished 2019)

Contains the following Python packges:
    
    qslam : Python package to implement NMQA algorithm. 
        Supports References [1], [2], [3].
 
    paduaq : Python package for Lagrange 2D interpolation at Padua points. 
        Supports Reference [2].
        Reference [3].

    clfanalysis : Python package for image classification for pre-processing
        single qubit measurements from trapped ions flouresence data.
        Supports Reference [1].

Contains the following directories for associated research analysis:

    expt_data : Converts trapped ion images into measurement database. 
        Codebase stored as a series of Jupyter Notebooks (.ipynb) 
        Supports Reference [1].
    
    expt_qslam : Compares NMQA/QSLAM to Naive Approach using experimental data.
        Codebase stored as a series of Python (.py) + Artemis (.pbs) scripts.
        Supports Reference [1].
    
    intrarun : Implements convergence analysis for NMQA-QSLAM.
        Codebase stored as a series of Python (.py) + Artemis (.pbs) scripts.
        Supports Reference [2].

    paperdata : Compares NMQA with Naive approaches using simulated data. 
        Codebase stored as a series of Python (.py) + Artemis (.pbs) scripts.
        Supports Reference [1].
   
    sampling_analysis: Compares NMQA with Padua interpolation techniques. 
        Codebase stored as a series of Python (.py) + Artemis (.pbs) scripts.
        Supports Reference [3].
    
    statedetectn : Implements image classification for pre-processing
        single qubit measurements from trapped ions flouresence data.
        Codebase stored as a series of Jupyter Notebooks (.ipynb).
        Supports Reference [1].
