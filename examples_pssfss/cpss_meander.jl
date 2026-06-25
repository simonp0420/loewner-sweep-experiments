include(joinpath(@__DIR__, "common.jl"))

function cpss_meander_case(; nfreq = 101, ntri = 600)
    outer(rot) = meander(;
        a = 3.97,
        b = 3.97,
        w1 = 0.13,
        w2 = 0.13,
        h = 2.53 + 0.13,
        units = mm,
        ntri,
        rot,
    )
    inner(rot) = meander(;
        a = 3.97 * sqrt(2),
        b = 3.97 / sqrt(2),
        w1 = 0.1,
        w2 = 0.1,
        h = 0.14 + 0.1,
        units = mm,
        ntri,
        rot,
    )
    center(rot) = meander(;
        a = 3.97,
        b = 3.97,
        w1 = 0.34,
        w2 = 0.34,
        h = 2.51 + 0.34,
        units = mm,
        ntri,
        rot,
    )

    substrate = Layer(width = 0.1mm, epsr = 2.6)
    foam(w) = Layer(width = w, epsr = 1.05)
    t1 = 4mm
    t2 = 2.45mm
    rot0 = 0

    strata = [
        Layer()
        outer(rot0)
        substrate
        foam(t1)
        inner(rot0 - 45)
        substrate
        foam(t2)
        center(rot0 - 2 * 45)
        substrate
        foam(t2)
        inner(rot0 - 3 * 45)
        substrate
        foam(t1)
        outer(rot0 - 4 * 45)
        substrate
        Layer()
    ]

    return strata, collect(range(10.0, 20.0, length = nfreq)), (θ = 0, ϕ = 0)
end

function run(; nfreq = 101, ntri = 600, seed = 1, run_direct = false, options = DEFAULT_SWEEP_OPTIONS)
    strata, freqs, steering = cpss_meander_case(; nfreq, ntri)
    return run_sweep_comparison(
        strata,
        freqs,
        steering;
        label = "Meanderline CPSS",
        options,
        seed,
        run_direct,
    )
end

results = run(run_direct = haskey(ENV, "RUN_DIRECT"));
plot_comparison(results; savepath = joinpath(@__DIR__, "cpss_meander.png"))
