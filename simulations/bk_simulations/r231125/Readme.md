# Goal

using selfing to model inbreeding


# Command

```
~/.local/bin/nextflow /local/chib/toconnor_grp/bing/posseleff_simulations/main.nf \
    -profile sge \
    -entry WF_MP \
    -resume \
    --mp_sets_json /local/chib/toconnor_grp/bing/posseleff_simulations/simulations/r231125/mp_genome_sets.json \
    --ibdne_no_diploid_convertion true \
    --peak_validate_meth ihs \
    --num_reps 1
    --resdir resdir_mp
```

# version:
```
commit 60f71cb4b8e579974e6c28cf880e4389d32ca60d (HEAD -> inbreeding2)
Author: Bing Guo <gbinux@gmail.com>
Date:   Sat Nov 25 11:16:44 2023 -0500

    allow selfing
```
