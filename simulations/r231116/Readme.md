# Goal

Try to rerun simulation similar to r230420 but with 
(1) ihs for peak validation
(2) no ibdne diploid conversion
(3) only remove high impact peaks for ibdne analysis to mitigate chromosome
fragmentation due to noise low-impact peaks.

# Simulation parameters

Please check `main.nf`. The variable `sp_sets` contains all parameters used for
single-population model, and `mp_sets` for parameters used for
multiple-population model.

# Command

```
nextflow /local/chib/toconnor_grp/bing/posseleff_simulations/main.nf \
    -profile sge \
    --num_reps 30 \
    --peak_validate_meth ihs \
    --ibdne_no_diploid_convertion true
```
