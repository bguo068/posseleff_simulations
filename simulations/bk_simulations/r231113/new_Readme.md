# Goal
Clean re-run r231113 for inbreeding effects

# Command

```
~/.local/bin/nextflow /local/chib/toconnor_grp/bing/posseleff_simulations/main.nf \
    -profile sge \
    -resume \
    --sp_sets_json ./genome_sets_sp.json \
    --mp_sets_json ./genome_sets_mp.json \
    --ibdne_no_diploid_convertion true \
    --peak_validate_meth ihs \
    --num_reps 1 
```

# verions: 

1.`posseleff_simulations`

```
commit 9c7c00f40b962e7d050a00c249e9645ddd9f95cf (HEAD -> inbreeding2)
Author: Bing Guo <gbinux@gmail.com>
Date:   Mon Nov 20 22:31:22 2023 -0500

    add notes on r231116 simulation
```

2. `ibdutils`
```
commit 35fac2212e112301c162926b3e1e540593713798 (HEAD -> main, origin/main, origin/HEAD)
Author: Bing Guo <gbinux@gmail.com>
Date:   Tue Nov 21 10:31:18 2023 -0500

    comment out dbg code
```
