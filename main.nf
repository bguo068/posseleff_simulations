#! /usr/bin/env nextflow

nextflow.enable.types = true

params {
    // Which simulations to run: 'all', 'sp' (single population) or 'mp'
    // (multiple population). Replaces the deprecated `-entry WF_SP` and
    // `-entry WF_MP` command line options.
    wf: String = 'all'

    ifm_transform: String = 'square' // 'square', 'cube' or 'none'
    ifm_ntrials: Integer = 1000
    test: Boolean = false
    num_reps: Integer = 1

    ibdne_mincm: Float = 2.0
    ibdne_minregion: Float = 10.0
    ibdne_flatmeth: String = 'none'
    min_ibdne_chr_min_afreq: Float = 0.2
    peak_validate_meth: String = 'ihs' // 'xirs' or 'ihs'
    sp_sets_json: String = ''
    mp_sets_json: String = ''

    // comma separated set numbers used to subseting
    // mp_sets/sp_sets when only some of them are needed.
    // example:  --genome_sets "10000,10003"
    genome_sets: String = ''

    // whether to skip ibdne and infomap analysis parts
    skip_ibdne: Boolean = false
    skip_ifm: Boolean = false

    // default to true so it is consistent with previous version
    // can be set to false on the nextflow command line
    ibdne_no_diploid_convertion: Boolean = true
    sim_multi_chrom_genome: Boolean = false
}

// A named set of simulation parameters for one genome (14 chromosomes)
record SimSet {
    label: String
    args: Map
}

// Ne estimates of a genome, before or after removing IBD peaks
record NeResult {
    label: String
    genome_set_id: Integer
    are_peaks_removed: Boolean
    ne: Path
}

// Infomap membership of a genome, before or after removing IBD peaks
record InfomapResult {
    label: String
    genome_set_id: Integer
    are_peaks_removed: Boolean
    member: Path
}

def sp_defaults() -> Map {
    return [
        seqlen : 100 * 15000,
        selpos : Math.round(0.33 * 100 * 15000),
        num_origins : 1,
        N : 10000,
        h : 0.5,
        s : 0.3,
        g_sel_start : 80,
        r : 0.01 / 15_000,
        sim_relatedness : 0,
        g_ne_change_start : 200,
        N0 : 1000,
        u : 1e-8,
        nsam : 1000, // haploid
    ]
}

def mp_defaults() -> Map {
    return [
        seqlen : 100 * 15000,
        selpos : Math.round(100 * 15000 * 0.33),
        num_origins : 1,
        N : 10000,
        h : 0.5,
        s : 0.3,
        g_sel_start : 80,
        r : 0.01 / 15_000,
        sim_relatedness : 0,
        mig : 1e-5,
        sel_mig : 0.01,
        npop : 5,
        nsam : 200, // haploid
        Tsplit : 500,
        u : 1e-8,
    ]
}

def sp_sets_builtin() -> Map {
    def d = sp_defaults()
    return [
        sp_neu: d + [s: 0.0, genome_set_id: 10000],
        sp_s01: d + [s: 0.1, genome_set_id: 10001],
        sp_s02: d + [s: 0.2, genome_set_id: 10002],
        sp_s03: d + [s: 0.3, genome_set_id: 10003],
        sp_g040:d + [g_sel_start: 40, genome_set_id: 10004],
        sp_g080:d + [g_sel_start: 80, genome_set_id: 10005],
        sp_g120:d + [g_sel_start:120, genome_set_id: 10006],
        sp_o01: d + [num_origins: 1, genome_set_id: 10007],
        sp_o03: d + [num_origins: 3, genome_set_id: 10008],
        sp_o27: d + [num_origins: 27, genome_set_id: 10009],
    ]
}

def mp_sets_builtin() -> Map {
    def d = mp_defaults()
    return [
        mp_s00: d + [s:0.0, genome_set_id: 20000],
        mp_s01: d + [s:0.1, genome_set_id: 20001],
        mp_s02: d + [s:0.2, genome_set_id: 20002],
        mp_s03: d + [s:0.3, genome_set_id: 20003],
    ]
}

// Render a map of arguments as command line options: --key1 value1 --key2 value2
def to_cli_args(args: Map) -> String {
    return args.keySet().collect { k -> "--${k} ${args[k]}" }.join(" ")
}

// Convert a map of genome sets (label -> simulation args) to a list of SimSet records
def to_sim_sets(sets: Map) -> List<SimSet> {
    return sets.keySet().collect { k -> record(label: k as String, args: sets[k] as Map) }.toList()
}

// if json file is provide, use params from json files
def genome_sets_from_json(default_: Map, json_fn: String) -> List<SimSet> {
    def config = new groovy.json.JsonSlurper().parse(file(json_fn)) as Map
    return config.keySet().collect { k -> record(label: k as String, args: default_ + (config[k] as Map)) }.toList()
}

// Build the simulation sets: builtin or json sets, optionally subset by
// genome set ids, replicated `num_reps` times, and reduced to a single small
// set in test mode.
def make_sim_sets(
    defaults: Map, builtin: Map, json_fn: String, genome_sets: String,
    num_reps: Integer, test: Boolean, test_args: Map
) -> List<SimSet> {
    def sets = json_fn != '' ? genome_sets_from_json(defaults, json_fn) : to_sim_sets(builtin)

    if (genome_sets != '') {
        def ids = genome_sets.tokenize(',').collect { id -> id.strip().toInteger() }.toList()
        sets = sets.findAll { s -> (s.args.genome_set_id as Integer) in ids }.toList()
    }

    if (num_reps != 1) {
        sets = (0..<num_reps).collectMany { rep ->
            sets.collect { s ->
                def gsid = (s.args.genome_set_id as Integer) + 10 * rep
                record(label: "${s.label}_rep${rep}", args: s.args + [genome_set_id: gsid])
            }
        }.toList()
    }

    if (test) {
        sets = sets.take(1).collect { s -> record(label: s.label, args: s.args + test_args) }.toList()
    }
    return sets
}

def chrno_from_file(f: Path) -> Integer {
    return f.getSimpleName().tokenize('_')[-1].toInteger()
}


process SIM_SP_CHR {
    tag "${args.genome_set_id}_${chrno}"

    input:
    record(label: String, chrno: Integer, args: Map)

    output:
    record(
        label: label,
        chrno: chrno,
        args: args,
        genome_set_id: args.genome_set_id as Integer,
        trees: file('*.trees'),
        vcf: file('*.vcf.gz'),
        daf: file('*.daf'),
        restart_count: file('*.restart_count'),
        true_ne: file('*.true_ne'),
    )

    script:
    def args_chr = to_cli_args(args + [chrno: chrno])
    """
    sim_single_pop.py $args_chr

    mkdir tmp; mv tmp_* tmp/
    """
    stub:
    def prefix = "${args.genome_set_id}_${chrno}"
    """
    touch ${prefix}{.trees,.vcf.gz,.daf,.restart_count,.true_ne}
    """
}

// for inbreeding simulation
process SIM_SP_GW {
    tag "${args.genome_set_id}_genome_wide"

    input:
    record(label: String, args: Map)

    output:
    record(
        label: label,
        args: args,
        genome_set_id: args.genome_set_id as Integer,
        trees: files('*.trees'),
        vcf: files('*.vcf.gz'),
        daf: file('*.daf'),
        restart_count: file('*.restart_count'),
        true_ne: file('*.true_ne'),
    )

    script:
    def args_chr = to_cli_args(args)
    """
    sim_single_pop_whole_genome.py $args_chr

    mkdir tmp; mv tmp_* tmp/
    """
    stub:
    def prefix = "${args.genome_set_id}"
    """
    touch ${prefix}{.daf,.restart_count,.true_ne}
    touch ${prefix}_{1..14}{.trees,.vcf.gz}
    """
}

process SIM_MP_CHR {
    tag "${args.genome_set_id}_${chrno}"

    input:
    record(label: String, chrno: Integer, args: Map)

    output:
    record(
        label: label,
        chrno: chrno,
        args: args,
        genome_set_id: args.genome_set_id as Integer,
        trees: file('*.trees'),
        vcf: file('*.vcf.gz'),
        daf: file('*.daf'),
        restart_count: file('*.restart_count'),
        demog: file('*_demog.png'),
    )

    script:
    def args_chr = to_cli_args(args + [chrno: chrno])
    """
    sim_multiple_pop.py $args_chr

    mkdir tmp; mv tmp_* tmp/
    """
    stub:
    def prefix = "${args.genome_set_id}_${chrno}"
    """
    touch ${prefix}{.trees,.vcf.gz,.daf,.restart_count,_demog.png}
    """
}

process CALL_IBD {
    tag "${args.genome_set_id}_${chrno}"

    input:
    record(label: String, chrno: Integer, args: Map, trees: Path, vcf: Path)

    output:
    record(
        label: label,
        chrno: chrno,
        args: args,
        genome_set_id: args.genome_set_id as Integer,
        vcf: vcf,
        ibd: file('*_tskibd.ibd'),
    )

    script:
    def args_local = to_cli_args([
        tree: trees,
        vcf: vcf,
        chrno: chrno,
        r: args.r,
        seqlen: args.seqlen,
        genome_set_id: args.genome_set_id,
        mincm: 2.0, call_hmmibd: 0, n:100, m:5,
    ])
    """
    call_ibd.py ${args_local}
    """
    stub:
    def prefix = "${args.genome_set_id}_${chrno}"
    """
    touch ${prefix}{.map,_tskibd.ibd}
    """
}

process PROC_DIST_NE {
    tag "${genome_set_id}"

    input:
    record(
        label: String,
        args: Map,
        genome_set_id: Integer,
        ibd_lst: List<Path>,
        vcf_lst: List<Path>,
        chrlist: String
    )
    ne_opts: Map

    output:
    record(
        label: label,
        args: args,
        genome_set_id: genome_set_id,
        ibdne_jar: file('ibdne.jar'),
        orig_sh: file('*_orig.sh'),
        orig_map: file('*_orig.map'),
        orig_ibd: file('*_orig.ibd.gz'),
        rmpeaks_sh: file('*_rmpeaks.sh'),
        rmpeaks_map: file('*_rmpeaks.map'),
        rmpeaks_ibd: file('*_rmpeaks.ibd.gz'),
        ibddist_ibd_obj: file('*.ibddist.ibdobj.gz'),
        ibdne_ibd_obj: files('*.ibdne.ibdobj.gz'),
    )

    script:
    def args_local = to_cli_args([
        ibd_files: ibd_lst.join(' '),
        vcf_files: vcf_lst.join(' '),
        genome_set_id: genome_set_id,
        ibdne_chrlist: chrlist,
    ] + ne_opts)

    """
    proc_dist_ne.py ${args_local}
    """
    stub:
    """
    touch ibdne.jar
    touch ${genome_set_id}{_orig.sh,_orig.map,_orig.ibd.gz}
    touch ${genome_set_id}{_rmpeaks.sh,_rmpeaks.map,_rmpeaks.ibd.gz}
    touch ${genome_set_id}.ibddist.ibdobj.gz
    touch ${genome_set_id}_orig.ibdne.ibdobj.gz
    touch ${genome_set_id}_rmpeaks.ibdne.ibdobj.gz
    """
}

process DAF_TO_CHRLIST {
    tag "${label}"

    input:
    record(label: String, daf: Path)
    min_afreq: Float

    output:
    record(label: label, chrlist: env('CHRLIST'))

    script:
    """
    CHRLIST=`daf2chrlist.py --daf $daf  --min_afreq $min_afreq`
    """
    stub:
    """
    CHRLIST="1 2 3"
    """
}

process PROC_INFOMAP {
    tag "${genome_set_id}"

    input:
    record(
        label: String,
        args: Map,
        genome_set_id: Integer,
        ibd_lst: List<Path>,
        vcf_lst: List<Path>
    )
    peak_validate_meth: String

    output:
    record(
        label: label,
        args: args,
        genome_set_id: genome_set_id,
        ifm_orig_ibd_obj: file('*_orig.ifm.ibdobj.gz'),
        ifm_rmpeaks_ibd_obj: file('*_rmpeaks.ifm.ibdobj.gz'),
    )

    script:
    def args_local = to_cli_args([
        ibd_files: ibd_lst.join(' '),
        vcf_files: vcf_lst.join(' '),
        genome_set_id: genome_set_id,
        peak_validate_meth: peak_validate_meth,
    ])
    """
    proc_infomap.py ${args_local}
    """
    stub:
    """
    touch ${genome_set_id}{_orig.ifm.ibdobj.gz,_rmpeaks.ifm.ibdobj.gz}
    """
}

process RUN_IBDNE {
    tag "${genome_set_id}_${are_peaks_removed}"

    input:
    record(
        label: String,
        genome_set_id: Integer,
        ibdne_jar: Path,
        ibdne_sh: Path,
        gmap: Path,
        ibd_gz: Path,
        are_peaks_removed: Boolean
    )

    output:
    record(
        label: label,
        genome_set_id: genome_set_id,
        are_peaks_removed: are_peaks_removed,
        ne: file('*.ne'),
    )

    script:
    """
    bash ${ibdne_sh}
    """
    stub:
    def src = are_peaks_removed ? "rmpeaks" : "orig"
    """
    touch ${genome_set_id}_${src}.ne
    """
}

process RUN_INFOMAP {
    tag "${args.genome_set_id}_${are_peaks_removed}"

    input:
    record(label: String, args: Map, ibd_obj: Path, are_peaks_removed: Boolean)
    ifm_opts: Map

    output:
    record(
        label: label,
        genome_set_id: args.genome_set_id as Integer,
        are_peaks_removed: are_peaks_removed,
        member: file('*_member.pq'),
    )

    script:
    def cut_mode = are_peaks_removed ? 'rmpeaks' : 'orig'
    def args_local = to_cli_args([
        ibd_obj: ibd_obj,
        npop: args.npop,
        nsam: args.nsam,
        genome_set_id: args.genome_set_id,
        cut_mode: cut_mode,
    ] + ifm_opts)
    """
    run_infomap.py ${args_local}
    """
    stub:
    def cut_mode = are_peaks_removed ? 'rmpeaks' : 'orig'
    """
    touch ${args.genome_set_id}_${cut_mode}_member.pq
    """
}



workflow WF_SP {
    take:
    sp_sets: Channel<SimSet>
    sim_multi_chrom_genome: Boolean
    skip_ibdne: Boolean
    min_chr_afreq: Float
    ne_opts: Map

    main:
    def nchrom = 14

    // SIMULATE chromosomes, either independently (one task per chromosome)
    // or as a whole genome (one task per genome) for inbreeding simulation
    ch_sets_chr = sp_sets.filter { _s -> !sim_multi_chrom_genome }
    ch_sets_gw = sp_sets.filter { _s -> sim_multi_chrom_genome }

    ch_sims_chr = SIM_SP_CHR(
        ch_sets_chr.flatMap { s ->
            (1..nchrom).collect { chrno -> record(label: s.label, chrno: chrno, args: s.args) }
        }
    )

    ch_sims_gw = SIM_SP_GW(ch_sets_gw)

    // split genome-wide outputs into per-chromosome trees/vcf, sorted by chromosome numbers
    ch_trees_vcf_gw = ch_sims_gw.flatMap { r ->
        r.trees.toSorted { t -> chrno_from_file(t) }.collect { t ->
            def chrno = chrno_from_file(t)
            def vcf = r.vcf.find { v -> chrno_from_file(v) == chrno }
            record(label: r.label, chrno: chrno, args: r.args, trees: t, vcf: vcf)
        }
    }
    ch_trees_vcf_chr = ch_sims_chr.map { r ->
        record(label: r.label, chrno: r.chrno, args: r.args, trees: r.trees, vcf: r.vcf)
    }

    // chr list of DAF > min_chr_afreq for genome-wide simulations, or all
    // chromosomes when chromosomes are simulated independently
    ch_chrlist = DAF_TO_CHRLIST(ch_sims_gw, min_chr_afreq)
        .mix(ch_sets_chr.map { s -> record(label: s.label, chrlist: (1..nchrom).join(" ")) })

    // CALL IBD segment per chromosome
    ch_ibd = CALL_IBD(ch_trees_vcf_chr.mix(ch_trees_vcf_gw))

    // COLLECT PER SETS
    ch_ibd_per_genome = ch_ibd
        .map { r -> tuple(r.label, nchrom, r) }
        .groupBy()
        .map { label, rs ->
            def sorted = rs.toSorted { r -> r.chrno }
            record(
                label: label,
                args: sorted[0].args,
                genome_set_id: sorted[0].genome_set_id,
                ibd_lst: sorted.collect { r -> r.ibd }.toList(),
                vcf_lst: sorted.collect { r -> r.vcf }.toList(),
            )
        }
        .join(ch_chrlist, by: 'label') // add chrlist
        .filter { _r -> !skip_ibdne }

    // Process IBD for ibd distribution and ne analyses
    ch_ne_input = PROC_DIST_NE(ch_ibd_per_genome, ne_opts)

    // RUN IBDNe actually
    // use mix so the second part do not have to wait for first part to finish complely
    ch_ne = RUN_IBDNE(
        ch_ne_input.map { r ->
            record(
                label: r.label, genome_set_id: r.genome_set_id, ibdne_jar: r.ibdne_jar,
                ibdne_sh: r.orig_sh, gmap: r.orig_map, ibd_gz: r.orig_ibd, are_peaks_removed: false,
            )
        }.mix(ch_ne_input.map { r ->
            record(
                label: r.label, genome_set_id: r.genome_set_id, ibdne_jar: r.ibdne_jar,
                ibdne_sh: r.rmpeaks_sh, gmap: r.rmpeaks_map, ibd_gz: r.rmpeaks_ibd, are_peaks_removed: true,
            )
        })
    )

    emit:
    sims: Channel<Record> = ch_sims_chr
    sims_gw: Channel<Record> = ch_sims_gw
    ibd: Channel<Record> = ch_ibd
    ne_input: Channel<Record> = ch_ne_input
    ne: Channel<NeResult> = ch_ne
}


workflow WF_MP {
    take:
    mp_sets: Channel<SimSet>
    skip_ifm: Boolean
    peak_validate_meth: String
    ifm_opts: Map

    main:
    def nchrom = 14

    // SIMULATE chromosomes
    ch_sims = SIM_MP_CHR(
        mp_sets.flatMap { s ->
            (1..nchrom).collect { chrno -> record(label: s.label, chrno: chrno, args: s.args) }
        }
    )

    // CALL IBD segment per chromosome
    ch_ibd = CALL_IBD(ch_sims)

    // COLLECT PER SETS
    ch_ibd_per_genome = ch_ibd
        .map { r -> tuple(r.label, nchrom, r) }
        .groupBy()
        .map { label, rs ->
            def sorted = rs.toSorted { r -> r.chrno }
            record(
                label: label,
                args: sorted[0].args,
                genome_set_id: sorted[0].genome_set_id,
                ibd_lst: sorted.collect { r -> r.ibd }.toList(),
                vcf_lst: sorted.collect { r -> r.vcf }.toList(),
            )
        }
        .filter { _r -> !skip_ifm }

    // Process IBD for population structure analyses
    ch_ifm_input = PROC_INFOMAP(ch_ibd_per_genome, peak_validate_meth)

    // RUN INFOMAP
    ch_ifm = RUN_INFOMAP(
        ch_ifm_input.map { r ->
            record(label: r.label, args: r.args, ibd_obj: r.ifm_orig_ibd_obj, are_peaks_removed: false)
        }.mix(ch_ifm_input.map { r ->
            record(label: r.label, args: r.args, ibd_obj: r.ifm_rmpeaks_ibd_obj, are_peaks_removed: true)
        }),
        ifm_opts
    )

    emit:
    sims: Channel<Record> = ch_sims
    ibd: Channel<Record> = ch_ibd
    ifm_input: Channel<Record> = ch_ifm_input
    ifm: Channel<InfomapResult> = ch_ifm
}

workflow {
    main:
    if (!(params.wf in ['all', 'sp', 'mp'])) {
        error "Invalid --wf '${params.wf}': expected 'all', 'sp' or 'mp'"
    }

    def sp_sets = params.wf in ['all', 'sp']
        ? make_sim_sets(
            sp_defaults(), sp_sets_builtin(), params.sp_sets_json, params.genome_sets,
            params.num_reps, params.test, [nsam: 50]
        )
        : []
    def mp_sets = params.wf in ['all', 'mp']
        ? make_sim_sets(
            mp_defaults(), mp_sets_builtin(), params.mp_sets_json, params.genome_sets,
            params.num_reps, params.test, [nsam: 50, npop: 2, N: 500]
        )
        : []

    def ne_opts = [
        peak_validate_meth: params.peak_validate_meth,
        ibdne_mincm: params.ibdne_mincm,
        ibdne_minregion: params.ibdne_minregion,
        ibdne_flatmeth: params.ibdne_flatmeth,
        ibdne_no_diploid_conversion: params.ibdne_no_diploid_convertion,
    ]
    def ifm_opts = [
        ntrials: params.ifm_ntrials,
        transform: params.ifm_transform,
    ]

    sp = WF_SP(
        channel.fromList(sp_sets),
        params.sim_multi_chrom_genome,
        params.skip_ibdne,
        params.min_ibdne_chr_min_afreq,
        ne_opts,
    )
    mp = WF_MP(
        channel.fromList(mp_sets),
        params.skip_ifm,
        params.peak_validate_meth,
        ifm_opts,
    )

    publish:
    sp_sims = sp.sims
    sp_gw_sims = sp.sims_gw
    mp_sims = mp.sims
    tskibd = sp.ibd.mix(mp.ibd)
    ne_input = sp.ne_input
    ne_output = sp.ne
    ifm_input = mp.ifm_input
    ifm_output = mp.ifm
}

// All outputs are published under `outputDir` (`results` by default, see
// nextflow.config), one subfolder per genome set: <genome_set_id>_<label>/
output {
    sp_sims {
        path { r ->
            r.trees >> "${r.genome_set_id}_${r.label}/trees/"
            r.vcf >> "${r.genome_set_id}_${r.label}/vcf/"
            r.daf >> "${r.genome_set_id}_${r.label}/daf/"
            r.restart_count >> "${r.genome_set_id}_${r.label}/restart_count/"
            r.true_ne >> "${r.genome_set_id}_${r.label}/true_ne/"
        }
    }

    sp_gw_sims {
        path { r ->
            r.trees >> "${r.genome_set_id}_${r.label}/trees/"
            r.vcf >> "${r.genome_set_id}_${r.label}/vcf/"
            r.daf >> "${r.genome_set_id}_${r.label}/daf/"
            r.restart_count >> "${r.genome_set_id}_${r.label}/restart_count/"
            r.true_ne >> "${r.genome_set_id}_${r.label}/true_ne/"
        }
    }

    mp_sims {
        path { r ->
            r.trees >> "${r.genome_set_id}_${r.label}/trees/"
            r.vcf >> "${r.genome_set_id}_${r.label}/vcf/"
            r.restart_count >> "${r.genome_set_id}_${r.label}/restart_count/"
            r.daf >> "${r.genome_set_id}_${r.label}/daf/"
            r.demog >> "${r.genome_set_id}_${r.label}/demog/"
        }
    }

    tskibd {
        path { r -> r.ibd >> "${r.genome_set_id}_${r.label}/tskibd/" }
    }

    ne_input {
        path { r ->
            [r.orig_sh, r.orig_map, r.orig_ibd, r.rmpeaks_sh, r.rmpeaks_map, r.rmpeaks_ibd] >> "${r.genome_set_id}_${r.label}/ne_input/"
            r.ibddist_ibd_obj >> "${r.genome_set_id}_${r.label}/ibddist_ibd/"
            r.ibdne_ibd_obj >> "${r.genome_set_id}_${r.label}/ibdne_ibd/"
        }
    }

    ne_output {
        path { r -> r.ne >> "${r.genome_set_id}_${r.label}/ne_output/" }
        index {
            path 'ne_output.csv'
            header true
        }
    }

    ifm_input {
        path { r ->
            [r.ifm_orig_ibd_obj, r.ifm_rmpeaks_ibd_obj] >> "${r.genome_set_id}_${r.label}/ifm_input/"
        }
    }

    ifm_output {
        path { r -> r.member >> "${r.genome_set_id}_${r.label}/ifm_output/" }
        index {
            path 'ifm_output.csv'
            header true
        }
    }
}
