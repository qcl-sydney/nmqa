# NMQA-QSLAM 

Research Codebase Noise Mapping for Quantum Architectures (NMQA). Formerly known as: Quantum Simultaneous Localisation and Mapping (QSLAM).

References:

[1] Swaroop Gupta, R., Milne, A. R., Edmunds, C. L., Hempel, C., & Biercuk, M. J. (2019). Adaptive scheduling of noise characterization in quantum computers. arXiv preprint arXiv:1904.07225. Accessed https://arxiv.org/abs/1904.07225 

[2] Gupta, R. S., & Biercuk, M. J. (2019). Convergence analysis for autonomous adaptive learning applied to quantum architectures. arXiv preprint arXiv:1911.05752. Accessed https://arxiv.org/abs/1911.05752

[3] Gupta, R., Govia L.~C.~G. & Biercuk, M. J. Interpolation and architectural impacts of spectator qubits for efficient quantum computer calibration and tuneup (forthcoming 2020).

Contains the following Python packges:
    
    qslam : Python package to implement NMQA-QSLAM algorithm. 
        Supports References [1], [2], [3].
 
    paduaq : Python package for Lagrange 2D interpolation at Padua points. 
        Supports Reference [3].

    clfanalysis : Python package for image classification for pre-processing
        single qubit measurements from trapped ion camera data.
        Supports Reference [1].

Contains the following directories for associated research analysis:

    nmqa1 : Compares NMQA-QSLAM to Naive Approach using simulated data. 
        Codebase stored as a series of Python (.py) + Artemis (.pbs) scripts.
        Supports Reference [1].
    
    nmqa2 : Implements convergence analysis for NMQA-QSLAM.
        Codebase stored as a series of Python (.py) + Artemis (.pbs) scripts.
        Supports Reference [2].
 
    nmqa3: Compares NMQA-QSLAM to Naive Approach with Padua interpolation. 
        Codebase stored as a series of Python (.py) + Artemis (.pbs) scripts.
        Supports Reference [3].
        
    expt_qslam : Compares NMQA-QSLAM to Naive Approach using experimental data.
        Codebase stored as a series of Python (.py) + Artemis (.pbs) scripts.
        Supports Reference [1].

    expt_data : Converts trapped ion images into measurement database. 
        Codebase stored as a series of Jupyter Notebooks (.ipynb) 
        Supports Reference [1].
    
    statedetectn : Analyses qubit state detection via image classification
        on trapped ion camera data.
        Codebase stored as a series of Jupyter Notebooks (.ipynb).
        Supports Reference [1].

