# Goal

simulate high inbreeding by assortative mating
for single pop population 



# Command

```
python3 generate_json.py

nextflow ../../main.nf \
    -profile sge  \
    -entry WF_SP \
    -resume \
    --sp_sets_json sp_genome_sets.json \
    --ibdne_no_diploid_convertion true \
    --peak_validate_meth ihs \
    --num_reps 1 
```

# Note: 

Non-default bypass parameters were used to allow higher inbreeding