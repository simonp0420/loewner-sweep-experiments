# Several parts were written directly using the PSSFSS internal code.
# I tried using PSSFSS as a black box and it worked,
# but I also wanted to test its performance, not just its accuracy.

using GLMakie
using LinearAlgebra
using Printf
using Random
using Statistics
using PSSFSS

include(joinpath(@__DIR__, "..", "src", "Loewner.jl"))
include(joinpath(@__DIR__, "..", "src", "random.jl"))
include(joinpath(@__DIR__, "..", "src", "randomLU.jl"))

const DEFAULT_SWEEP_OPTIONS = SweepOptions(
    Δf = 1e-5,
    q1 = 8,
    q2 = 12,
    adaptive = false,
    tol = 1e-3,
    l = 0.04,
    p = 4,
    memory = 3,
    use_D = false,
    data_partition = true,
    parallel = true,
    nthreads = Threads.nthreads() ÷ 2
)

function gsm_matrix(result)
    gsm = result.gsm
    return Matrix{ComplexF64}([gsm[1, 1] gsm[1, 2]; gsm[2, 1] gsm[2, 2]])
end

db20(x) = 20 * log10(max(abs(x), eps(Float64)))
sparam_db(mats, i, j) = [db20(mat[i, j]) for mat in mats]

function relative_errors(approx, reference)
    return [norm(approx[i] - reference[i]) / max(norm(reference[i]), eps(Float64)) for i in eachindex(reference)]
end

function pssfss_evaluated_freqs(logfile)
    isfile(logfile) || return Float64[]
    freqs = Float64[]
    for line in eachline(logfile)
        m = match(r"^\s*([0-9]+(?:\.[0-9]*)?(?:[eE][+-]?[0-9]+)?) GHz\s*$", line)
        m === nothing && continue
        push!(freqs, parse(Float64, m.captures[1]))
    end
    return sort!(unique(freqs))
end

as_float_vector(x) = x isa Number ? Float64[float(x)] : Float64.(float.(collect(x)))

#PSSFSS function
function pssfss_internal_setup(strata, freqs, steering)
    layers = PSSFSS.Layer[deepcopy(s) for s in strata if s isa PSSFSS.Layer]
    sheets = PSSFSS.RWGSheet[s for s in strata if s isa PSSFSS.Sheet]

    islayer = map(x -> x isa PSSFSS.Layer, strata)
    issheet = map(x -> x isa PSSFSS.RWGSheet, strata)
    nl = length(layers)
    nj = nl - 1
    ns = length(sheets)

    sint = cumsum(islayer)[issheet]
    junc = zeros(Int, nj)
    junc[sint] = 1:ns

    stkeys = Tuple(keys(steering))
    stvalues = [as_float_vector(s) for s in steering]
    all(length.(stvalues) .== 1) || error("Internal adapter supports one steering point only")
    steer = PSSFSS.getsttuple(stkeys, stvalues[1][1], stvalues[2][1])

    PSSFSS.check_inputs(layers, sheets, junc, freqs, stkeys, stvalues, [])

    k0min, k0max = PSSFSS.twopi * 1e9 / PSSFSS.c₀ .* extrema(freqs)
    gbls = PSSFSS.choose_gblocks(layers, sheets, junc, k0min)
    gsm_save = Vector{PSSFSS.GSM}(undef, length(gbls))
    PSSFSS.choose_layer_modes!(layers, sheets, junc, gbls, k0max, PSSFSS.dbmin)
    gbldup = PSSFSS.get_gbldup(gbls, layers, sheets, junc)
    usi = PSSFSS.unique_indices(sheets)
    rwgdat = PSSFSS.RWGData[PSSFSS.setup_rwg(sheet) for sheet in sheets]

    uvec = map(sheets) do sh
        sh.style == "NULL" && return 0.0
        ufactor = 0.5 * ustrip(Float64, sh.units, 1u"m")
        ufactor * max(norm(sh.β₁), norm(sh.β₂))
    end

    β00k1 = if keys(steer)[1] == :ψ₁
        ψ1, ψ2 = (deg2rad(x) for x in steer)
        upm = ustrip(Float64, sheets[1].units, 1u"m")
        β1, β2 = sheets[1].β₁ * upm, sheets[1].β₂ * upm
        (ψ1 * β1 + ψ2 * β2) / PSSFSS.twopi
    else
        θ, ϕ = steer
        st = sind(θ)
        sp, cp = sincosd(ϕ)
        PSSFSS.SVector(st * cp, st * sp)
    end

    return (; layers, sheets, junc, freqs, steer, β00k1, gbls, gsm_save, gbldup, usi, rwgdat, uvec)
end


function make_pssfss_internal_solver(strata, freqs, steering)
    setup = pssfss_internal_setup(strata, freqs, steering)

    solve(fGHz) = begin
        out = PSSFSS.compute_next_freq(
            fGHz,
            setup.β00k1,
            setup.steer,
            setup.layers,
            setup.sheets,
            setup.usi,
            setup.rwgdat,
            setup.uvec,
            setup.junc,
            setup.gbls,
            setup.gbldup,
            setup.gsm_save,
        )
        out === nothing && error("PSSFSS skipped $(fGHz) GHz due to cutoff principal modes")
        _, result = out
        gsm_matrix(result)
    end

    return solve, setup
end

function run_sweep_comparison(strata, freqs, steering; label, options = DEFAULT_SWEEP_OPTIONS, seed = 1, run_random = true, run_direct = false,)

    seed === nothing || Random.seed!(seed)

    println("Example:            ", label)
    println("Frequency grid:     ", length(freqs), " points from ", first(freqs), " to ", last(freqs), " GHz")


    setup_time = @elapsed begin
        solve, setup = make_pssfss_internal_solver(strata, freqs, steering)
    end

    println("Building Loewner model and evaluating surrogate...")

    tloewner_build = @elapsed model = sweep(solve, freqs, options)
    tloewner_eval = @elapsed approx = model.(freqs)
    tloewner = setup_time + tloewner_build + tloewner_eval

    if run_random
        println("Building random Loewner model and evaluating surrogate...")
        trandom_build = @elapsed random_model = sweep_random(solve, freqs, options)
        trandom_eval = @elapsed random_approx = random_model.(freqs)
        trandom = setup_time + trandom_build + trandom_eval
    end

    println("Running PSSFSS default fastsweep at all grid points...")
    fast_log = tempname() * "_pssfss_fast.log"
    tfast = @elapsed fast_results = analyze(
        strata,
        freqs,
        steering;
        showprogress = false,
        logfile = fast_log,
        resultfile = devnull,
    )
    fast = gsm_matrix.(fast_results)
    fast_f0 = pssfss_evaluated_freqs(fast_log)

    direct = nothing
    tdirect = NaN
    relerr = nothing
    relerr_random = nothing
    relerr_fast = nothing
    if run_direct
        println("Running direct reference at all grid points...")
        tdirect = @elapsed direct = solve.(freqs)
        relerr = relative_errors(approx, direct)
        relerr_random = run_random ? relative_errors(random_approx, direct) : nothing
        relerr_fast = relative_errors(fast, direct)
    end

    println()
    @printf("Internal setup time: %.3f s\n", setup_time)
    @printf("Loewner total time:  %.3f s\n", tloewner)
    @printf("  build time:        %.3f s\n", tloewner_build)
    @printf("  eval time:         %.3f s\n", tloewner_eval)
    if run_random
        @printf("Random total time:   %.3f s\n", trandom)
        @printf("  build time:        %.3f s\n", trandom_build)
        @printf("  eval time:         %.3f s\n", trandom_eval)
    end
    run_direct && @printf("Direct grid time:    %.3f s\n", tdirect)
    @printf("PSSFSS fast time:    %.3f s\n", tfast)
    println("Loewner samples:     ", length(model.f₀), " / ", length(freqs))
    run_random && println("Random samples:      ", length(random_model.f₀), " / ", length(freqs))
    println("PSSFSS fast samples: ", length(fast_f0), " / ", length(freqs))
    println("PSSFSS fast log:     ", fast_log)

    if run_direct
        imax = argmax(relerr)
        imax_fast = argmax(relerr_fast)
        @printf("Loewner max relerr:  %.3e at %.3f GHz\n", relerr[imax], freqs[imax])
        @printf("Loewner med relerr:  %.3e\n", median(relerr))
        if run_random
            imax_random = argmax(relerr_random)
            @printf("Random max relerr:   %.3e at %.3f GHz\n", relerr_random[imax_random], freqs[imax_random])
            @printf("Random med relerr:   %.3e\n", median(relerr_random))
        end
        @printf("PSSFSS fast maxerr:  %.3e at %.3f GHz\n", relerr_fast[imax_fast], freqs[imax_fast])
        @printf("PSSFSS fast mederr:  %.3e\n", median(relerr_fast))
    end

    return (;
        label,
        freqs,
        steering,
        direct,
        approx,
        random_approx,
        fast,
        relerr,
        relerr_random,
        relerr_fast,
        model,
        random_model,
        fast_f0,
        fast_log,
        setup,
        setup_time,
        tloewner,
        tloewner_build,
        tloewner_eval,
        trandom,
        trandom_build,
        trandom_eval,
        tdirect,
        tfast,
    )
end

# Just works
function plot_comparison(results; i = 1, j = 1, dense_surrogate = false, dense_n = 1001, savepath = joinpath(@__DIR__, "comparison.png"), legend_outside = true, legend_position = :rt,)

    fplot = dense_surrogate ?
        collect(range(first(results.freqs), last(results.freqs), length = dense_n)) :
        results.freqs
    surrogate_plot = dense_surrogate ? results.model.(fplot) : results.approx
    random_plot = dense_surrogate && results.random_model !== nothing ?
        results.random_model.(fplot) :
        results.random_approx

    direct_db = results.direct === nothing ? nothing : sparam_db(results.direct, i, j)
    surrogate_db = sparam_db(surrogate_plot, i, j)
    random_db = random_plot === nothing ? nothing : sparam_db(random_plot, i, j)
    fast_db = sparam_db(results.fast, i, j)
    sample_db = sparam_db(results.model.H₀, i, j)
    random_sample_db = results.random_model === nothing ? nothing : sparam_db(results.random_model.H₀, i, j)

    fast_sample_idx = [findfirst(==(f), results.freqs) for f in results.fast_f0]
    fast_sample_idx = filter(!isnothing, fast_sample_idx)
    fast_sample_freqs = results.freqs[fast_sample_idx]
    fast_sample_db = fast_db[fast_sample_idx]

    fig = Figure(size = legend_outside ? (1120, 620) : (1000, 620))
    ax = Axis(
        fig[1, 1];
        xlabel = "Frequency (GHz)",
        ylabel = "Magnitude (dB)",
        title = "$(results.label): S$(i)$(j)",
    )

    lines!(ax, results.freqs, fast_db; linewidth = 2.4, color = :black, linestyle = :solid, label = "PSSFSS fastsweep")
    direct_db !== nothing &&
        lines!(ax, results.freqs, direct_db; linewidth = 2.5, color = :gray35, label = "PSSFSS direct")
    lines!(
        ax,
        fplot,
        surrogate_db;
        linewidth = 2.2,
        color = :dodgerblue3,
        linestyle = :dash,
        label = dense_surrogate ? "Loewner dense" : "Loewner",
    )
    random_db !== nothing &&
        lines!(ax, fplot, random_db; linewidth = 2.2, color = :darkorange2, linestyle = :dashdot, label = dense_surrogate ? "Random dense" : "Random")
    scatter!(
        ax,
        results.model.f₀,
        sample_db;
        markersize = 12,
        color = :red,
        strokewidth = 0.8,
        strokecolor = :white,
        label = "Loewner samples",
    )
    if random_sample_db !== nothing
        scatter!(
            ax,
            results.random_model.f₀,
            random_sample_db;
            markersize = 13,
            color = :magenta4,
            marker = :utriangle,
            strokewidth = 0.8,
            strokecolor = :white,
            label = "Random samples",
        )
    end
    scatter!(
        ax,
        fast_sample_freqs,
        fast_sample_db;
        markersize = 9,
        color = :seagreen3,
        marker = :diamond,
        strokewidth = 0.8,
        strokecolor = :white,
        label = "PSSFSS fast samples",
    )

    if legend_outside
        Legend(fig[1, 2], ax; framevisible = true)
        colgap!(fig.layout, 12)
    else
        axislegend(ax; position = legend_position, framevisible = true)
    end
    if savepath !== nothing
        save(savepath, fig)
        println("Saved figure to: ", savepath)
    end
    display(fig)
    return fig
end
