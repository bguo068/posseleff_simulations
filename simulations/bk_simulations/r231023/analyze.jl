using PyCall
using CairoMakie
using DataFrames
using CSV
using Parquet2
using Clustering # has randindex function to calc indices including adjusted rand index
using CategoricalArrays # conver string to categorical
using NamedArrays
using Statistics

# use svg, so vscode plot pane has better preview

CairoMakie.activate!(type="svg", pt_per_unit=0.75)


struct SimParams
    label::String
    gsid::Int
    s::Float64
    rel::Int
    power::Float64
    delta::Float64
    relg::Int
    SimParams(
        label::String,
        gsid::Int,
        s::Float64,
        rel::Int,
        power::Float64,
        delta::Float64,
        relg::Int,
    ) = new(label, gsid, s, rel, power, delta, relg)
end

# constructor with defaults
SimParams(label::String, gsid::Int, s::Float64) = SimParams(label, gsid, s, 0, 1.0, 0.01, 40)



const SIMS = [
    # SimParams("sp_neu", 10000, 0.0),
    # SimParams("sp_s03", 10003, 0.3),
    SimParams("sp_rels03_11a", 111000, 0.3, 1, 10.0, 10.0, 200),
    SimParams("sp_rels00_11a", 110001, 0.0, 1, 10.0, 10.0, 200),
    SimParams("sp_rels03_11b", 110002, 0.3, 1, 10.0, 10.0, 100),
    SimParams("sp_rels00_11b", 110003, 0.0, 1, 10.0, 10.0, 100),
    SimParams("sp_rels03_11c", 110004, 0.3, 1, 10.0, 10.0, 50),
    SimParams("sp_rels00_11c", 110005, 0.0, 1, 10.0, 10.0, 50),
    SimParams("sp_rels03_12a", 120000, 0.3, 1, 10.0, 1.0, 200),
    SimParams("sp_rels00_12a", 120001, 0.0, 1, 10.0, 1.0, 200),
    SimParams("sp_rels03_12b", 120002, 0.3, 1, 10.0, 1.0, 100),
    SimParams("sp_rels00_12b", 120003, 0.0, 1, 10.0, 1.0, 100),
    SimParams("sp_rels03_12c", 120004, 0.3, 1, 10.0, 1.0, 50),
    SimParams("sp_rels00_12c", 120005, 0.0, 1, 10.0, 1.0, 50),
    SimParams("sp_rels03_13a", 130000, 0.3, 1, 10.0, 0.01, 200),
    SimParams("sp_rels00_13a", 130001, 0.0, 1, 10.0, 0.01, 200),
    SimParams("sp_rels03_13b", 130002, 0.3, 1, 10.0, 0.01, 100),
    SimParams("sp_rels00_13b", 130003, 0.0, 1, 10.0, 0.01, 100),
    SimParams("sp_rels03_13c", 130004, 0.3, 1, 10.0, 0.01, 50),
    SimParams("sp_rels00_13c", 130005, 0.0, 1, 10.0, 0.01, 50),
    SimParams("sp_rels03_14a", 140000, 0.3, 1, 10.0, 1e-8, 200),
    SimParams("sp_rels00_14a", 140001, 0.0, 1, 10.0, 1e-8, 200),
    SimParams("sp_rels03_14b", 140002, 0.3, 1, 10.0, 1e-8, 100),
    SimParams("sp_rels00_14b", 140003, 0.0, 1, 10.0, 1e-8, 100),
    SimParams("sp_rels03_14c", 140004, 0.3, 1, 10.0, 1e-8, 50),
    SimParams("sp_rels00_14c", 140005, 0.0, 1, 10.0, 1e-8, 50),
    SimParams("sp_rels03_21a", 211000, 0.3, 1, 1.0, 10.0, 200),
    SimParams("sp_rels00_21a", 210001, 0.0, 1, 1.0, 10.0, 200),
    SimParams("sp_rels03_21b", 210002, 0.3, 1, 1.0, 10.0, 100),
    SimParams("sp_rels00_21b", 210003, 0.0, 1, 1.0, 10.0, 100),
    SimParams("sp_rels03_21c", 210004, 0.3, 1, 1.0, 10.0, 50),
    SimParams("sp_rels00_21c", 210005, 0.0, 1, 1.0, 10.0, 50),
    SimParams("sp_rels03_22a", 220000, 0.3, 1, 1.0, 1.0, 200),
    SimParams("sp_rels00_22a", 220001, 0.0, 1, 1.0, 1.0, 200),
    SimParams("sp_rels03_22b", 220002, 0.3, 1, 1.0, 1.0, 100),
    SimParams("sp_rels00_22b", 220003, 0.0, 1, 1.0, 1.0, 100),
    SimParams("sp_rels03_22c", 220004, 0.3, 1, 1.0, 1.0, 50),
    SimParams("sp_rels00_22c", 220005, 0.0, 1, 1.0, 1.0, 50),
    SimParams("sp_rels03_23a", 230000, 0.3, 1, 1.0, 0.01, 200),
    SimParams("sp_rels00_23a", 230001, 0.0, 1, 1.0, 0.01, 200),
    SimParams("sp_rels03_23b", 230002, 0.3, 1, 1.0, 0.01, 100),
    SimParams("sp_rels00_23b", 230003, 0.0, 1, 1.0, 0.01, 100),
    SimParams("sp_rels03_23c", 230004, 0.3, 1, 1.0, 0.01, 50),
    SimParams("sp_rels00_23c", 230005, 0.0, 1, 1.0, 0.01, 50),
    SimParams("sp_rels03_24a", 240000, 0.3, 1, 1.0, 1e-8, 200),
    SimParams("sp_rels00_24a", 240001, 0.0, 1, 1.0, 1e-8, 200),
    SimParams("sp_rels03_24b", 240002, 0.3, 1, 1.0, 1e-8, 100),
    SimParams("sp_rels00_24b", 240003, 0.0, 1, 1.0, 1e-8, 100),
    SimParams("sp_rels03_24c", 240004, 0.3, 1, 1.0, 1e-8, 50),
    SimParams("sp_rels00_24c", 240005, 0.0, 1, 1.0, 1e-8, 50),
    SimParams("sp_rels03_31a", 311000, 0.3, 1, 0.1, 10.0, 200),
    SimParams("sp_rels00_31a", 310001, 0.0, 1, 0.1, 10.0, 200),
    SimParams("sp_rels03_31b", 310002, 0.3, 1, 0.1, 10.0, 100),
    SimParams("sp_rels00_31b", 310003, 0.0, 1, 0.1, 10.0, 100),
    SimParams("sp_rels03_31c", 310004, 0.3, 1, 0.1, 10.0, 50),
    SimParams("sp_rels00_31c", 310005, 0.0, 1, 0.1, 10.0, 50),
    SimParams("sp_rels03_32a", 320000, 0.3, 1, 0.1, 1.0, 200),
    SimParams("sp_rels00_32a", 320001, 0.0, 1, 0.1, 1.0, 200),
    SimParams("sp_rels03_32b", 320002, 0.3, 1, 0.1, 1.0, 100),
    SimParams("sp_rels00_32b", 320003, 0.0, 1, 0.1, 1.0, 100),
    SimParams("sp_rels03_32c", 320004, 0.3, 1, 0.1, 1.0, 50),
    SimParams("sp_rels00_32c", 320005, 0.0, 1, 0.1, 1.0, 50),
    SimParams("sp_rels03_33a", 330000, 0.3, 1, 0.1, 0.01, 200),
    SimParams("sp_rels00_33a", 330001, 0.0, 1, 0.1, 0.01, 200),
    SimParams("sp_rels03_33b", 330002, 0.3, 1, 0.1, 0.01, 100),
    SimParams("sp_rels00_33b", 330003, 0.0, 1, 0.1, 0.01, 100),
    SimParams("sp_rels03_33c", 330004, 0.3, 1, 0.1, 0.01, 50),
    SimParams("sp_rels00_33c", 330005, 0.0, 1, 0.1, 0.01, 50),
    SimParams("sp_rels03_34a", 340000, 0.3, 1, 0.1, 1e-8, 200),
    SimParams("sp_rels00_34a", 340001, 0.0, 1, 0.1, 1e-8, 200),
    SimParams("sp_rels03_34b", 340002, 0.3, 1, 0.1, 1e-8, 100),
    SimParams("sp_rels00_34b", 340003, 0.0, 1, 0.1, 1e-8, 100),
    SimParams("sp_rels03_34c", 340004, 0.3, 1, 0.1, 1e-8, 50),
    SimParams("sp_rels00_34c", 340005, 0.0, 1, 0.1, 1e-8, 50),
]

# const RESDIR::String = "/local/chib/toconnor_grp/bing/posseleff_simulations/simulations/r231019/resdir/"
const RESDIR::String = "resdir/"


# prepare data --------------------
function get_ibd_obj_fname(resdir::String, gsid::Int, label::String)::String
    p::String = "$(resdir)/$(gsid)_$(label)/ibddist_ibd/$(gsid)_2.0_10.0_none.ibddist.ibdobj.gz"
    @assert isfile(p) "path $(p) does not exist"
    return p
end

function get_ibd_obj_fname2(resdir::String, gsid::Int, label::String)::String
    p::String = "$(resdir)/$(gsid)_$(label)/ibdne_ibd/$(gsid)_2.0_10.0_none_orig.ibdne.ibdobj.gz"
    @assert isfile(p) "path $(p) does not exist"
    return p
end

function get_ne_file(resdir::String, gsid::Int, label::String, rm::String)::String
    p::String = "$(resdir)/$(gsid)_$(label)/ne_output/$(gsid)_2.0_10.0_none_$(rm).ne"
    @assert isfile(p) "path $(p) does not exist"
    return p
end

function get_true_ne_file(resdir::String, gsid::Int, label::String)::String
    p::String = "$(resdir)/$(gsid)_$(label)/true_ne/$(gsid).true_ne"
    @assert isfile(p) "path $(p) does not exist"
    return p
end

function read_ne(ne_file)::DataFrame
    df = DataFrame(CSV.File(ne_file))
    rename!(df, Symbol.(["GEN", "NE", "L95", "U95"]))
    df
end

function read_true_ne(ne_file)::DataFrame
    df = DataFrame(CSV.File(ne_file)) # with two columns ["GEN", "NE"]
    df
end

function get_daf(resdir::String, gsid::Int, label::String)::String
    p::String = "$(resdir)/$(gsid)_$label/daf/$(gsid).daf"
    @assert isfile(p) "path $p does not exist"
    return p
end

function read_daf(fname::String)::DataFrame
    DataFrame(CSV.File(fname))
end

# python function
py"""
from ibdutils.utils.ibdutils import IBD
def read_ibd_obj_py(ibdobjfname):
    ibd = IBD.pickle_load(ibdobjfname)
    return (
        ibd._df,
        ibd._cov_df,
        ibd._peaks_df,
        ibd._genome._chr_df,
        ibd._genome._drg_df,
        ibd._genome._gmap.gmap,
    )
def get_unrelated_sample_count(ibdobjfname):
    ibd = IBD.pickle_load(ibdobjfname)
    mat = ibd.make_ibd_matrix()
    unrelated_samples = ibd.get_unrelated_samples(mat)
    return unrelated_samples.shape[0]
"""
"""helper function convert pandas dataframe (PyObject) to 
julia DataFrame """
function pd_to_df(df_pd::Union{PyObject,Nothing})
    df = DataFrame()
    if df_pd === nothing
        return df
    end
    for col in df_pd.columns
        df[!, col] = getproperty(df_pd, col).values
    end
    df
end

const IbdObj = NamedTuple{(:ibd, :cov, :peaks, :chrs, :drg, :gmap),NTuple{6,DataFrame}}

# julia function
function read_ibd_obj(ibdobjfname::String)::IbdObj
    # get a 6 tuple
    tup6 = py"read_ibd_obj_py"(ibdobjfname)
    # convert to julia DataFrame
    return (
        ibd=pd_to_df(tup6[1]),
        cov=pd_to_df(tup6[2]),
        peaks=pd_to_df(tup6[3]),
        chrs=pd_to_df(tup6[4]),
        drg=pd_to_df(tup6[5]),
        gmap=pd_to_df(tup6[6]),
    )
end
function read_ibd_obj_cov_only(ibdobjfname::String)
    # get a 6 tuple
    tup6 = py"read_ibd_obj_py"(ibdobjfname)
    # convert to julia DataFrame
    return (pd_to_df(tup6[2]))
end
function read_ibd_obj_ibd_only(ibdobjfname::String)::DataFrame
    # get a 6 tuple
    tup6 = py"read_ibd_obj_py"(ibdobjfname)
    # convert to julia DataFrame
    return pd_to_df(tup6[1])
end
function read_ibd_obj_n_unrelated(ibdobjfname::String)::Int
    return py"get_unrelated_sample_count"(ibdobjfname)
end

function get_total_mat(ibd::DataFrame)::Matrix{Float64}
    gdf = groupby(ibd, [:Id1, :Id2])
    tot = combine(gdf, [:Start, :End] => ((s, e) -> sum((e .- s) ./ 15000.0)) => :TotIbd)
    n::Integer = max(maximum(tot.Id1), maximum(tot.Id2)) + 1
    m = zeros((n, n))
    for r in eachrow(tot)
        m[r.Id1+1, r.Id2+1] = r.TotIbd
        m[r.Id2+1, r.Id1+1] = r.TotIbd
    end
    return m
end

function clust_total_mat(m::Matrix{Float64})::Vector{Int}
    mi = minimum(m)
    ma = maximum(m)
    m2 = ma ./ (m .- mi .+ 1e-3)
    o = hclust(m2, :average).order
    return o
end


function create_axes_and_add_labels()
    f = Figure(resolution=(1200 * 2, 800 * 2))
    axes = Dict{NTuple{4,Int},Axis}()
    for i in 1:4
        for j in 1:3
            for inner_i in 1:2
                for inner_j in 1:3
                    ax = Axis(f[i, j][inner_i, inner_j])
                    ax.xticksvisible = false
                    ax.yticksvisible = false
                    ax.xticklabelsvisible = false
                    ax.yticklabelsvisible = false
                    if inner_i == 1
                        g = [200, 100, 50][inner_j]
                        ax.title = "g=$g"
                    end
                    if inner_j == 1
                        s = [0.0, 0.3][inner_i]
                        ax.ylabel = "s=$s"
                    end
                    axes[(i, j, inner_i, inner_j)] = ax
                end
            end
        end
    end
    for i in 1:4
        Box(f[i, 0], color=:gray90)
        delta = [10.0, 1.0, 0.01, 1e-8][i]
        label = "delta=$(delta)"
        Label(f[i, 0], label, rotation=pi / 2, tellheight=false)
    end
    for j in 1:3
        Box(f[0, j], color=:gray90)
        power = [10.0, 1.0, 0.1][j]
        label = "power=$(power)"
        Label(f[0, j], label, tellwidth=false)
    end
    return (f, axes)
end

# calculate, cluster and plot total IBD map
function get_ax(axes::Dict{NTuple{4,Int},Axis}, sim::SimParams)::Axis
    i = findfirst(x -> x == sim.delta, [10.0, 1.0, 0.01, 1e-8])
    j = findfirst(x -> x == sim.power, [10.0, 1.0, 0.1])
    inner_i = findfirst(x -> x == sim.s, [0.0, 0.3])
    inner_j = findfirst(x -> x == sim.relg, [200, 100, 50])
    return axes[(i, j, inner_i, inner_j)]
end


function plot_totibd(f::Figure, axes::Dict{NTuple{4,Int},Axis})
    for sim in SIMS
        println(sim.gsid, ' ', sim.label)
        ax = get_ax(axes, sim)
        try
            fname = get_ibd_obj_fname(RESDIR, sim.gsid, sim.label)
            ibd = read_ibd_obj_ibd_only(fname)
            m = get_total_mat(ibd)
            o = clust_total_mat(m)
            m2 = m[o, o]
            heatmap!(ax, m2, colormap=:Blues)
            n_unrel = read_ibd_obj_n_unrelated(fname)
            text!(ax, 0, 0, text="num_unrel = $n_unrel", align=(:left, :top))
        catch e
            println(e)
        end
    end
    save("out_totibd_mat_clust.png", f)
end

function plot_ne(f::Figure, axes::Dict{NTuple{4,Int},Axis})
    for sim::SimParams in SIMS
        println(sim.gsid, ' ', sim.label)
        ax = get_ax(axes, sim)
        # plot true ne
        try
            ne_fname = get_true_ne_file(RESDIR, sim.gsid, sim.label)
            df = read_true_ne(ne_fname)
            lines!(ax, df.GEN, df.NE, label="Truth", color=:black, linestyle=:dot)
            # 
            for rm in ["orig", "rmpeaks"]
                try
                    ne_fname = get_ne_file(RESDIR, sim.gsid, sim.label, rm)
                    df = read_ne(ne_fname)
                    color = :blue
                    ls = :solid
                    if rm == "rmpeaks"
                        color = :red
                        ls = :dot
                    end
                    # here using a factor 4 to account for the fact that we are using the
                    # haploid genomes as homozygous diploids 
                    lines!(ax, df.GEN, df.NE / 4.0, color=color, label=rm, linestyle=ls)
                catch e
                    println(e)
                end
            end
        catch e
            println(e)
        end
        xlims!(ax, 0, 100)
        ylims!(ax, 1e2, 1e6)
        ax.yscale = log10
        # axislegend()
    end
    save("out_ne.png", f)
end

function plot_cov(f::Figure, axes::Dict{NTuple{4,Int},Axis})
    for sim in SIMS
        println(sim.gsid, ' ', sim.label)
        ax = get_ax(axes, sim)
        try
            cov = read_ibd_obj_cov_only(get_ibd_obj_fname(RESDIR, sim.gsid, sim.label))
            println("mean", sim, mean(cov.Coverage))
            lines!(ax, cov.GwStart, cov.Coverage)
            # ylims!(ax, 0, 1e5)
            linkyaxes!(ax, axes[(1, 1, 1, 1)])
        catch e
            println(e)
        end
    end
    save("out_cov.png", f)
end

function plot_daf(f::Figure, axes::Dict{NTuple{4,Int},Axis})
    for sim in SIMS
        println(sim.gsid, ' ', sim.label)
        ax = get_ax(axes, sim)
        try
            daf = read_daf(get_daf(RESDIR, sim.gsid, sim.label))
            num_established_selection = 0
            for i in 1:14
                col = "DAF_CHR$i"
                if maximum(daf[!, col]) > 0.2
                    num_established_selection += 1
                end
                lines!(ax, daf.GEN, daf[!, col])

                linkyaxes!(ax, axes[(1, 1, 1, 1)])
            end
            num_peaks = 0
            try
                num_peaks = size(read_ibd_obj(get_ibd_obj_fname2(RESDIR, sim.gsid, sim.label)).peaks, 1)
            catch e
                println(e)
            end
            text!(ax, 0, 0, text="est=$(num_established_selection), peaks=$(num_peaks)", align=(:left, :bottom))
        catch e
            println(e)
        end
    end
    save("out_daf.png", f)
end

f, axes = create_axes_and_add_labels()
plot_cov(f, axes)

f, axes = create_axes_and_add_labels()
plot_totibd(f, axes)

f, axes = create_axes_and_add_labels()
plot_ne(f, axes)

#
f, axes = create_axes_and_add_labels()
plot_daf(f, axes)
#

# scratch space
"""
tmp = read_ibd_obj_ibd_only(get_ibd_obj_fname2(RESDIR, SIMS[end-1].gsid, SIMS[end-1].label))
unique(tmp, [:Id1, :Id2])
unique(hcat(tmp.Id1, tmp.Id2))
"""
