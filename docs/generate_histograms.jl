using ESRI
using Plots
using Random
using SparseArrays

const OUTDIR = joinpath(@__DIR__, "src", "assets")

default(
    fmt = :svg,
    size = (760, 420),
    legend = false,
    grid = :y,
    linewidth = 0,
)

function build_demo_economy()
    Random.seed!(42)
    n = 1_000
    W = sprand(n, n, 0.01)
    W[1:n+1:end] .= 0
    info = IndustryInfo(rand(1:4, n), [true, true, false, false])
    return ESRIEconomy(W, info), W
end

function plot_scores(values, title, path; color = :steelblue)
    bins = length(values) < 20 ? max(length(values), 5) : 16
    p = histogram(
        values;
        bins = bins,
        color = color,
        alpha = 0.85,
        xlabel = "ESRI score",
        ylabel = "Number of firms",
        title = title,
        linecolor = :white,
        framestyle = :box,
        left_margin = 6Plots.mm,
        bottom_margin = 5Plots.mm,
    )
    savefig(p, path)
end

function plot_score_comparison(values_a, values_b, label_a, label_b, title, path)
    p = histogram(
        values_a;
        bins = 16,
        normalize = :pdf,
        alpha = 0.45,
        color = :steelblue,
        label = label_a,
        xlabel = "ESRI score",
        ylabel = "Density",
        title = title,
        linecolor = :steelblue,
        framestyle = :box,
        left_margin = 6Plots.mm,
        bottom_margin = 5Plots.mm,
    )
    histogram!(
        p,
        values_b;
        bins = 16,
        normalize = :pdf,
        alpha = 0.45,
        color = :firebrick,
        label = label_b,
        linecolor = :firebrick,
    )
    savefig(p, path)
end

mkpath(OUTDIR)

econ, W = build_demo_economy()

scores = esri(econ; maxiter = 40, tol = 1e-3, threads = false)
plot_scores(scores, "Distribution of default firm-shock ESRI scores", joinpath(OUTDIR, "scores_hist.svg"))

output_weights = vec(sum(W; dims = 2))
scores_output = esri(
    econ;
    final_weights = output_weights,
    combine = :min,
    maxiter = 25,
    tol = 1e-3,
    threads = false,
)
spike_weights = copy(output_weights)
spike_weights[10] = 100.0
scores_spike = esri(
    econ;
    final_weights = spike_weights,
    combine = :min,
    maxiter = 25,
    tol = 1e-3,
    threads = false,
)
plot_score_comparison(
    scores_output,
    scores_spike,
    "Precomputed output weights",
    "Firm 10 weight set to 100",
    "Effect of changing one firm's final weight",
    joinpath(OUTDIR, "scores_custom_weights_compare.svg"),
)

subset_indices = collect(25:25:1_000)
subset_scores = esri(econ; firm_indices = subset_indices, maxiter = 20, tol = 1e-3, threads = false)
selected_scores = subset_scores[subset_scores .> 0]
plot_scores(
    selected_scores,
    "Distribution for a selected subset of 40 shocked firms",
    joinpath(OUTDIR, "subset_scores_hist.svg");
    color = :darkorange,
)
