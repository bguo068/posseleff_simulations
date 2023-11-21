# Goal

Try to rerun simulation similar to r230420 but with ihs for peak validation

# Change 

To keep the mp_sets and sp_sets consistent with those in r230420 as
they were changed since r230420, (commit: `334008cd55c9d40c76fbea0318c4991e0e8d1fe1`)

# Command

```
nextflow ../../main.nf \
    -profile sge \
    --num_reps 30 \
    --peak_validate_meth ihs \
    --ibdne_no_diploid_convertion true
```

# versions:

ibdutils @ 9ceffdccb8e856b256cf91ba21525c61328ce57b
posseleff_simulations @ db9a30bcb8ff65955eb4ebbe820a6d272523d2a9
