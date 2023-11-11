# Goal

- to run more repeats on top of r231106 (with selected genome sets to reduce number of simulations)
- to understand how randomness affect the overall relatedness 
- explore how selection start time affect the detection of selection effect

# command

```
~/.local/bin/nextflow /local/chib/toconnor_grp/bing/posseleff_simulations/main.nf \
    -entry WF_SP \
    -resume \
    -profile sge \
    --ibdne_no_diploid_convertion true \
    --num_reps 20 \
    --sim_inbreeding true

```

# findings

1. The effect of selection and bia-correction are not stable across 20 repeats
2. The averge trend is still concordant with expectation:
    - without inbreeding, selection-bias-correction is easily seen
    - with inbreeding, correction seems not working or over correcting.
3. selection time window of 50 generation might be good for inbreeding
simulation, but can be too short for non-inbreeding simulation
    - non-inbreeding and inbreeding cases might have different speed of 
      elevation of the allele frequencies of selected mutations
    - using different selection time length for noninbreeding and inbreeding
      simulation
