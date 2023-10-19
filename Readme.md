# `posseleff_simulations`

The `posseleff_simulations` repository is designed to assess the impact of
positive selection on Identity by Descent (IBD)-based inferences, leveraging
population genetic simulation and true IBD methodologies. The pipeline begins by running
simulations to generate tree sequences, allowing for both single and multiple
population models for various analyses. The selection simulation can be tailored
with parameters such as the selection coefficient, number of origins of the
favored mutation, selection starting time, and high-relatedness simulation
options. From there, the pipeline involves calling true IBD segments from the
tree sequence, followed by IBD processing for generating input files, selection
correction, and calling IBDNe and Infomap for Ne and population structure
inference. The pipeline is highly configurable, accommodating different
scenarios and requirements, and produces detailed output for further analysis.
Installation and execution instructions are provided, as well as options for
customization, making it a versatile tool for many IBD-related analyses.


## Outline of the pipeline
- Run simulation and generate tree sequence
    - Two demographic models
        - Single population model for IBD distribution and Ne analysis, and 
        - Multiple population model for population structure analysis
    - Selection simulation:
        - Each of the above models allows for positive selection simulation
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
    - Remove IBD within validated IBD peak region and generate a
    selection-corrected version of IBD for calling IBDNe and Infomap
- Call `IBDNe` and `Infomap` for Ne and population structure inference

## System requirements and Software environment

The pipeline has been tested on Linux Operation system and can be easily adapted to MacOS with
simple changes. Software dependencies and the version numbers are specified in the
'./env.yaml' Conda recipe. Additional depencies that are not available from Conda are
specified in the installation instruction below. The overall installation time
is about 5-15 minutes.

To create the software environment:
1. Install nextflow. See [nextflow documentation](https://www.nextflow.io/docs/latest/getstarted.html)
2. Install conda from [here](https://docs.conda.io/en/latest/miniconda.html) if you have not
3. Install software: `python3 ./init.py`, this will
    - Detect and Install `simulation` Conda environment
        - including [`ibdutils`](https://github.com/bguo068/ibdutils)
    - Download and compile [`tskibd`](https://github.com/bguo068/tskibd) for
    true IBD inference from tree sequence
    - Download
    [`IBDNe`](https://faculty.washington.edu/browning/ibdne/ibdne.23Apr20.ae9.jar)


## How to run the pipeline

1. Activate the `simulation` environment: `conda activate simulation`
2. Run the pipeline: `nextflow ./main.nf -profile sge --num_reps 30 -resume`. 
3. For large datasets, using a cluster such as SGE is recommended. An example
`sge` profile is provided in the `nextflow.config` file and should be adjusted
to fit your cluster system. If run on a local computer, please remove the
`-profile sge` option from the above command.
4. The pipeline can be reconfigured in the following files
    - Pipeline parameters can be found top lines in `main.nf`
    - More simulation parameters can be found in the definition of
     `sp_defaults`,`mp_defaults`, `sp_sets` and `mp_sets` dictionaries within
     the `main.nf` file

## Input, Test Data and Output

1. No input or test data is needed. If desired, pipeline can be reconfigured as
mentioned above. The pipeline can be tested by commenting all but one entry in
the `sp_sets` or `mp_sets`.
2. Output files/folders:
    - each subfolder of `resdir` represents simulations for a set of chromosomes
    - within each subfolder:
        - `ifm_output` contains infomap results
        - `ne_output` contains IBDNe estimates

# Citations

If you find this repository useful, please cite our preprint:
> Guo, B., Borda, V., Laboulaye, R., Spring, M. D., Wojnarski, M., Vesely, B.
A., Silva, J. C., Waters, N. C., O'Connor, T. D., & Takala-Harrison, S. (2023).
Strong Positive Selection Biases Identity-By-Descent-Based Inferences of Recent
Demography and Population Structure in Plasmodium falciparum. bioRxiv : the
preprint server for biology, 2023.07.14.549114.
https://doi.org/10.1101/2023.07.14.549114

Other citations:

- `Xir,s` statistics: 
> Henden, L., Lee, S., Mueller, I., Barry, A., & Bahlo, M. (2018).
Identity-by-descent analyses for measuring population dynamics and selection in
recombining pathogens. PLoS genetics, 14(5), e1007279.
https://doi.org/10.1371/journal.pgen.1007279

- `IBDNe`
> Browning, S. R., & Browning, B. L. (2015). Accurate Non-parametric Estimation
of Recent Effective Population Size from Segments of Identity by Descent.
American journal of human genetics, 97(3), 404–418.
https://doi.org/10.1016/j.ajhg.2015.07.012

- `Infomap` algorithm
> Rosvall, M., & Bergstrom, C. T. (2008). Maps of random walks on complex
networks reveal community structure. Proceedings of the National Academy of
Sciences of the United States of America, 105(4), 1118–1123.
https://doi.org/10.1073/pnas.0706851105


## Related Repository:
`tskibd`: https://github.com/bguo068/tskibd


## Inbreeding

Branch 1: `inbreeding` focusing on modifying  simulation structure for multiple chromosome genome simulation for (Single-population model)
branch 2: `inbreeding2` focusing on fine-tunning parameters for allow established selection
