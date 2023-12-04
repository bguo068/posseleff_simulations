# Goal

- compare iHS vs Xirs based validation of IBD peaks (due to positive selection)
- rerun pipeline similar to r231112 but with iHS for peak validation

# Command for simulation

1. share the same SIM_SP_CHR/CALL_IBD process with r231112 but differ in other processes
in which `ihs` method is used to validate peaks
```
~/.local/bin/nextflow /local/chib/toconnor_grp/bing/posseleff_simulations/main.nf \
    -entry WF_SP \
    -resume \
    -profile sge \
    --ibdne_no_diploid_convertion true \
    --peak_validate_meth ihs \
    --resdir resdir_sp \
    --num_reps 1 
```

2. share the same SIM_MP_CHR/CALL_IBD process with r231111 but differ in other processes
in which `ihs` method is used to validate peaks
```
~/.local/bin/nextflow /local/chib/toconnor_grp/bing/posseleff_simulations/main.nf \
    -profile sge \
    -entry WF_MP \
    -resume \
    --ibdne_no_diploid_convertion true \
    --peak_validate_meth ihs \
    --resdir resdir_mp \
    --num_reps 1 
```