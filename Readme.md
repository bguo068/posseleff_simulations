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
    - High-relatedness/inbreeding simulations
        - support inbreeding modeling via three strategies:
            * shrinking population size
            * positive assortative mating
            * selfing
        - see `simulations/Readme.md` for examples.
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
1. Install package manager `pixi` from [here](https://pixi.prefix.dev/latest/installation/) if you have not
2. Clone the repo: `git clone https://github.com/bguo068/posseleff_simulations.git`
3. Change the cloned directory: `cd posseleff_simulations`
4. Install packages according the `pixi` recipe by running `pixi install` 
5. Install additional software: `pixi run setup` (runs `./init.py`), this will
    - Download and compile [`tskibd`](https://github.com/bguo068/tskibd) for
    true IBD inference from tree sequence
    - Download
    [`IBDNe`](https://faculty.washington.edu/browning/ibdne/ibdne.23Apr20.ae9.jar)


## How to run the pipeline

1. Run the pipeline locally: `pixi run nextflow ./main.nf -resume`. Or, you can
2. Run the pipeline on hpc: `pixi run nextflow ./main.nf -profile sge --num_reps 30 -resume` or use `-profile slurm` for using slurm scheduler.
3. For large datasets, using a cluster such as SGE/Slurm is recommended. An example
`sge` profile is provided in the `nextflow.config` file and should be adjusted
to fit your cluster system. If run on a local computer, please remove the
`-profile sge` option from the above command.
4. The pipeline can be reconfigured in the following files
    - Pipeline parameters are declared in the typed `params` block at the top of `main.nf`
    - More simulation parameters can be found in the `sp_defaults()`,
     `mp_defaults()`, `sp_sets_builtin()` and `mp_sets_builtin()` functions
     within the `main.nf` file
    - For more complicated or large scale simulations, the `--sp_sets_json`,
    `--mp_sets_json` and `--wf` arguments are recommended.  see
    `simulations/Readme.md` for examples.
    - `--wf sp` or `--wf mp` runs only the single- or multiple-population
    simulations (default: `--wf all`). This replaces the `-entry WF_SP` /
    `-entry WF_MP` options used with older versions of the pipeline.
5. The pipeline requires Nextflow >= 26.10 (installed by `pixi`). It uses
static typing (`nextflow.enable.types`) and workflow outputs (the `output`
block in `main.nf`). Run `pixi run nextflow lint main.nf` to type-check changes.

## Input, Test Data and Output

1. No input or test data is needed. If desired, pipeline can be reconfigured as
mentioned above. The pipeline can be tested with `--test`, which runs only the
first entry of `sp_sets` and `mp_sets` with a small sample size, or quickly
dry-run with `-stub`.
2. Output files/folders (published as symlinks into the `outputDir` set in
`nextflow.config`, default `results`; override with `-output-dir <dir>`):
    - each subfolder of `results` represents simulations for a set of chromosomes
    - within each subfolder:
        - `ifm_output` contains infomap results
        - `ne_output` contains IBDNe estimates
    - `results/ne_output.csv` and `results/ifm_output.csv` are index files listing
    the IBDNe and infomap result files with their label, genome set id and
    whether IBD peaks were removed

# Citations

If you find this repository useful, please cite our preprint:
> Guo B, Borda V, Laboulaye R, Spring MD, Wojnarski M, Vesely BA, Silva JC, Waters NC, O'Connor TD, Takala-Harrison S. Strong positive selection biases identity-by-descent-based inferences of recent demography and population structure in Plasmodium falciparum. Nat Commun. 2024 Mar 20;15(1):2499. doi: 10.1038/s41467-024-46659-0. PMID: 38509066; PMCID: PMC10954658.

Other citations:

- `iHS` statistics and calculation:

iHS calculation via scikit-allel: 
> Miles, A. et al. cggh/scikit-allel: v1.3.7. (2023) doi:10.5281/ZENODO.8326460. 

iHS statistics:
> Voight, B. F., Kudaravalli, S., Wen, X. & Pritchard, J. K. A map of recent
positive selection in the human genome. PLoS Biol. 4, e72 (2006).

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
- `tskibd`: https://github.com/bguo068/tskibd
- `ibdutils`: https://github.com/bguo068/ibdutils
- empirical analysis pipeline: https://github.com/bguo068/posseleff_empirical


