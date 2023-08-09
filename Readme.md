# `posseleff_simulations`

## Outline of the pipeline
- Run simulation and generate tree sequence
    - Two demographic models
        - Single population model for IBD distribution and Ne analysis, and 
        - Multiple population model for population structure analysis
    - Selection simulation:
        - Each of the above models allows positive selection simulation
        - Tunable Simulation parameters include
            1.	selection coefficient
            2.	number of origins of the favored mutation
            3.	selection starting time (generations)
    - High-relatedness simulation (0: off, 1: on)
- Call true IBD segment from tree sequence
- IBD processing and selection correction
    - IBD processing for generating input files for IBDNe (Ne estimation)
    - IBD processing for generating IBD for calling Infomap (population structure)
    - Identify and validate IBD peaks (due to selection)
    - Remove IBD within validated IBD peak region and generate a selection-corrected version of IBD for calling IBDNe and Infomap
- Call `IBDNe` and `Infomap` for Ne and population structure inference


## How to run the pipeline

1. Install nextflow. See [nextflow documentation](https://www.nextflow.io/docs/latest/getstarted.html)
2. Install software: `python3 ./init.py`, this will
    - Install `simulation` Conda environment
        - including [`ibdutils`](https://github.com/bguo068/ibdutils)
    - Download and compile [`tskibd`](https://github.com/bguo068/tskibd) for
    true IBD inference from tree sequence
    - Download
    [`IBDNe`](https://faculty.washington.edu/browning/ibdne/ibdne.23Apr20.ae9.jar)
3. Activate the `simulation` environment: `conda activate simulation`
4. Run the pipeline: `nextflow ../main.nf -profile sge --num_reps 30 -resume`. 
5. Pipeline can be reconfigured in the following files
    - SGE cluster and resource allocation can be found `nextflow.config` file
    - Pipeline parameters can found top lines in `main.nf`
    - More simulation parameters can be found in the defintion of
     `sp_defaults`,`mp_defaults`, `sp_sets` and `mp_sets` dictionaries within `main.nf`

## Input and Output

1. No input needed. If desired, pipeline can be reconfigure as mentioned above.
2. Output files/folders:
    - each subfolder of `resdir` represents simulations for a set of chromosomes
    - within each subfolder:
        - `ifm_output` contains infomap results
        - `ne_output` contains IBDNe estimates
