#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

params.ifm_transform = ["square", "cube", "none"][0]
params.ifm_ntrails = 1000
params.test = false

def sp_defaults = [
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
    N0 : 3000,
    u : 1e-8,
    nsam : 100,
]

def mp_defaults = [
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
    nsam : 100,
    Tsplit : 500,
    u : 1e-8,
]

def sp_sets = [
    sp_neu: sp_defaults + [s: 0.0, genome_set_id: 10000],
    sp_s01: sp_defaults + [s: 0.1, genome_set_id: 10001],
    sp_s02: sp_defaults + [s: 0.2, genome_set_id: 10002],
    sp_s03: sp_defaults + [s: 0.3, genome_set_id: 10003],
    sp_g040:sp_defaults + [g_sel_start: 40, genome_set_id: 10004],
    sp_g080:sp_defaults + [g_sel_start: 80, genome_set_id: 10005],
    sp_g120:sp_defaults + [g_sel_start:120, genome_set_id: 10006],
    sp_o01: sp_defaults + [num_origins: 1, genome_set_id: 10007],
    sp_o03: sp_defaults + [num_origins: 3, genome_set_id: 10008],
    sp_o27: sp_defaults + [num_origins: 27, genome_set_id: 10009],
]

def mp_sets = [
    mp_s00: mp_defaults + [s:0.0, genome_set_id: 20000],
    mp_s01: mp_defaults + [s:0.1, genome_set_id: 20001],
    mp_s02: mp_defaults + [s:0.2, genome_set_id: 20002],
    mp_s03: mp_defaults + [s:0.3, genome_set_id: 20003],
]


process SIM_SP_CHR {
    input: 
    tuple val(label), val(chrno), val(args)

    output: 
    tuple val(label), val(chrno), path("*.trees"), path("*.vcf.gz"), emit: trees_vcf
    tuple val(label), val(chrno), path("*.daf"), emit: daf
    tuple val(label), val(chrno), path("*.restart_count"), emit: restart_count
    tuple val(label), val(chrno), path("*.true_ne"), emit: true_ne

    script:
    def args_chr = (args + [chrno: chrno]).collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    sim_single_pop.py $args_chr

    mkdir tmp; mv tmp_* tmp/
    """
    stub:
    def prefix="${args.genome_set_id}_${chrno}"
    """
    touch ${prefix}{.trees,.vcf.gz,.daf,.restart_count,.true_ne}
    """
}

process SIM_MP_CHR {
    input:
    tuple val(label), val(chrno), val(args)

    output:
    tuple val(label), val(chrno), path("*.trees"), path("*.vcf.gz"), emit: trees_vcf
    tuple val(label), val(chrno), path("*.daf"), emit: daf
    tuple val(label), val(chrno), path("*.restart_count"), emit: restart_count
    tuple val(label), val(chrno), path("*_demog.png"), emit: demog

    script:
    def args_chr = (args + [chrno: chrno]).collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    sim_multiple_pop.py $args_chr

    mkdir tmp; mv tmp_* tmp/
    """
    stub:
    def prefix="${args.genome_set_id}_${chrno}"
    """
    touch ${prefix}{.trees,.vcf.gz,.daf,.restart_count,_demog.png}
    """
}

process CALL_IBD {
    input: 
        tuple val(label), val(chrno), val(args),  path(trees), path(vcf)
    output:
        tuple val(label), val(chrno), path("*_tskibd.ibd"), emit: tskibd
    script:
    def args_local = [
        tree: trees,
        vcf: vcf,
        chrno: chrno,
        r: args.r,
        seqlen: args.seqlen,
        genome_set_id: args.genome_set_id,
        mincm: 2.0, call_hmmibd: 0, n:100, m:5,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    call_ibd.py ${args_local}
    """
    stub:
    def prefix="${args.genome_set_id}_${chrno}"
    """
    touch ${prefix}{.map,_tskibd.ibd}
    """
}

process PROC_DIST_NE {
    input:
        tuple val(label), path(ibd_lst), val(genome_set_id)
    output:
        tuple val(label), path("ibdne.jar"), path("*_orig.sh"), \
                path("*_orig.map"), path("*_orig.ibd.gz"), emit: ne_input_orig
        tuple val(label), path("ibdne.jar"), path("*_rmpeaks.sh"),  \
                path("*_rmpeaks.map"), path("*_rmpeaks.ibd.gz"), emit: ne_input_rmpeaks
        tuple val(label), path("*_ibddist_ibd.pq"), emit: ibddist_ibd_pq
    script:
    def args_local = [
        ibd_files: "${ibd_lst}", // path is a blank separate list
        genome_set_id: genome_set_id,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    proc_dist_ne.py ${args_local} 
    """
    stub:
    """
    touch ibdne.jar
    touch ${genome_set_id}{_orig.sh,_orig.map,_orig.ibd.gz}
    touch ${genome_set_id}{_rmpeaks.sh,_rmpeaks.map,_rmpeaks.ibd.gz}
    touch ${genome_set_id}_ibddist_ibd.pq
    """
}

process PROC_INFOMAP {
    input:
        tuple val(label), path(ibd_lst), val(genome_set_id)
    output:
        tuple val(label), path("*_ifm_orig_ibd.pq"), emit: ifm_orig_ibd_pq
        tuple val(label), path("*_ifm_rmpeaks_ibd.pq"), emit: ifm_rmpeaks_ibd_pq
    script:
    def args_local = [
        ibd_files: "${ibd_lst}", // path is a blank separate list
        genome_set_id: genome_set_id,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    proc_infomap.py ${args_local}
    """
    stub:
    """
    touch ${genome_set_id}{_ifm_orig_ibd.pq,_ifm_rmpeaks_ibd.pq}
    """
}

process RUN_IBDNE {
    input:
        tuple val(label), path(ibdne_jar), path(ibdne_sh), path(gmap), path(ibd_gz), \
            val(are_peaks_removed), val(args)
    output:
        tuple val(label), val(are_peaks_removed), path("*.ne")
    script:
    """
    bash ${ibdne_sh}
    """
    stub:
    def src = are_peaks_removed ? "rmpeaks": "orig"
    """
    touch ${args}_${src}.ne
    """
}

process RUN_INFOMAP {
    input:
        tuple val(label), path(ibd_pq), val(are_peaks_removed), val(args)
    output:
        tuple val(label), val(are_peaks_removed), path("*_member.pq")
    script:
    def args_local = [
        ibd_pq: ibd_pq,
        npop: args.npop,
        nsam: args.nsam,
        genome_set_id: args.genome_set_id,
        ntrails: params.ifm_ntrails,
        transform: params.ifm_transform,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    run_infomap.py ${args_local}
    """
    stub:
    """
    touch ${args.genome_set_id}_member.pq
    """
}



workflow WF_SP {
    // SIMULATION SETS

    println("\n\n------------- Single Population Sets-------------")
    sp_sets.each{k, v -> println("${k}\t\t${v}")}


    // SIMULATE chromosomes
    chr_chrno = channel.from(1..14)
    ch_sp_params = channel.from(sp_sets.collect{k, v-> [k, v]})

    if (params.test) {
        ch_sp_params = ch_sp_params.first().map{
            label, args -> def args2 = args + [nsam:50]; [label, args2]
        }
    }

    SIM_SP_CHR(ch_sp_params.combine(chr_chrno).map{a,b,c->[a,c,b]})
    // SIM_SP_CHR.out.trees_vcf.view{it-> 
    //    "${it[0]}\t${it[1]}\t${it[2].getName()}\t${it[3].getName()}"}



    // CALL IBD segment per chromosome
    CALL_IBD(SIM_SP_CHR.out.trees_vcf.combine(ch_sp_params, by:0)
        .map{label,chrno,trees,vcf,args-> [label,chrno, args, trees, vcf]}
    )
    // CALL_IBD.out.tskibd.view{label, chrno, ibd -> [label, chrno, ibd.getName()]}



    // COLLECT PER SETS
    ch_ibd_per_genome = CALL_IBD.out.tskibd
        .map{ label, chrno, ibd -> [groupKey(label, 14), [chrno, ibd]] }
        .groupTuple(by:0, sort: {x,y-> x[0]<=>y[0] })
        .combine(ch_sp_params, by:0)
        .map{label, ll, args-> [label, ll.collect{it[1]}, args.genome_set_id]}

    // ch_ibd_per_genome.view{label, trees_lst, genome_set_id->
    //        [label, trees_lst.collect{it.getName()}, genome_set_id]}


    // Process IBD for ibd distribution and ne analyses
    PROC_DIST_NE(ch_ibd_per_genome)
    // PROC_DIST_NE.out.ibddist_ibd_pq.view{label, ibdpq -> [label, ibdpq.getName()]}


    // RUN IBDNe actually
    RUN_IBDNE(
        PROC_DIST_NE.out.ne_input_orig.map{it -> it + false}.combine(ch_sp_params, by:0).concat(
            PROC_DIST_NE.out.ne_input_rmpeaks.map{it -> it + true}.combine(ch_sp_params, by:0)
        )
    )

    // RUN_IBDNE.out.map{label, yn, ne->[label, yn, ne.getName()]}.view()
    

}


workflow WF_MP {
    // SIMULATION SETS
    
    println("\n\n------------- Multiple Population Sets-------------")
    mp_sets.each{k, v -> println("${k}\t\t${v}")}



    // SIMULATE chromosomes
    chr_chrno = channel.from(1..14)
    ch_mp_params = channel.from(mp_sets.collect{k, v-> [k, v]})

    if (params.test) {
        ch_mp_params = ch_mp_params.first().map{
            label, args -> def args2 = args + [nsam:50, npop:2, N:500]; [label, args2]
        }
    }

    SIM_MP_CHR(ch_mp_params.combine(chr_chrno).map{a,b,c->[a,c,b]})
    // SIM_MP_CHR.out.trees_vcf.view{it-> 
    //    "${it[0]}\t${it[1]}\t${it[2].getName()}\t${it[3].getName()}"}



    // CALL IBD segment per chromosome
    CALL_IBD(SIM_MP_CHR.out.trees_vcf.combine(ch_mp_params, by:0)
        .map{label,chrno,trees,vcf,args-> [label,chrno, args, trees, vcf]}
    )
    // CALL_IBD.out.tskibd.view{label, chrno, ibd -> [label, chrno, ibd.getName()]}



    // COLLECT PER SETS
    ch_ibd_per_genome = CALL_IBD.out.tskibd
        .map{ label, chrno, ibd -> [groupKey(label, 14), [chrno, ibd]] }
        .groupTuple(by:0, sort: {x,y-> x[0]<=>y[0] })
        .combine(ch_mp_params, by:0)
        .map{label, ll, args-> [label, ll.collect{it[1]}, args.genome_set_id]}

    // ch_ibd_per_genome.view{label, trees_lst, genome_set_id->
    //        [label, trees_lst.collect{it.getName()}, genome_set_id]}


    // Process IBD for ibd distribution and ne analyses
    PROC_INFOMAP(ch_ibd_per_genome)
    
    // PROC_INFOMAP.out.ifm_orig_ibd_pq.view{label, ibdpq -> [label, ibdpq.getName()]}


    // RUN INFOMAP
    RUN_INFOMAP(
        PROC_INFOMAP.out.ifm_orig_ibd_pq.map{it -> it + false}.combine(ch_mp_params, by:0).concat(
            PROC_INFOMAP.out.ifm_rmpeaks_ibd_pq.map{it -> it + true}.combine(ch_mp_params, by:0)
        )
    )

    // RUN_INFOMAP.out.view{label, yns, member -> [label, yns, member.getName()]}

}

workflow {
    WF_SP()
    WF_MP()
}