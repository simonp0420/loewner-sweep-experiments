include(joinpath(@__DIR__, "common.jl"))

function splitring_case(; nfreq = 651, ntri = 450)
    r123 = [3.7, 4.25, 4.8]
    d = 10.8
    w = 0.3
    g = 0.3
    a = r123 .- w / 2
    b = a .+ w

    sheet = splitring(;
        class = 'M',
        units = mm,
        sides = 42,
        ntri,
        a,
        b,
        s1 = [d, 0],
        s2 = [0, d],
        gapwidth = [g, g, g],
        gapcenter = [-90, 90, -90],
    )

    strata = [
        Layer()
        sheet
        Layer(ϵᵣ = 10.2, tanδ = 0.0023, width = 0.13mm)
        Layer()
    ]

    return strata, collect(range(1.0, 14.0, length = nfreq)), (θ = 0, ϕ = 0)
end

function run(; nfreq = 651, ntri = 450, seed = 1, run_direct = false, options = DEFAULT_SWEEP_OPTIONS)
    strata, freqs, steering = splitring_case(; nfreq, ntri)
    return run_sweep_comparison(
        strata,
        freqs,
        steering;
        label = "Split-ring resonator",
        options,
        seed,
        run_direct,
    )
end

results = run();
plot_comparison(results; savepath = joinpath(@__DIR__, "splitring.png"))
