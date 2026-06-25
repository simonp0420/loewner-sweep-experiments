include(joinpath(@__DIR__, "common.jl"))

function jerusalem_wang_werner_case(; nfreq = 381, ntri = 800)
    P = 2.5
    L1 = 3P / 4
    L2 = P / 8
    w = L2
    A = P / 2
    B = P / 16

    sheet = jerusalemcross(;
        class = 'M',
        w,
        L1,
        L2,
        P,
        A,
        B,
        ntri,
        units = cm,
    )

    strata = [
        Layer()
        sheet
        Layer(epsr = 2.0, tandel = 0.05, width = 0.02cm)
        Layer()
    ]

    freqs = collect(range(1.0, 20.0, length = nfreq))
    steering = (phi = 1, theta = 45)
    return strata, freqs, steering
end

function run(; nfreq = 381, ntri = 800, seed = 1, run_direct = false, options = DEFAULT_SWEEP_OPTIONS)
    strata, freqs, steering = jerusalem_wang_werner_case(; nfreq, ntri)
    return run_sweep_comparison(
        strata,
        freqs,
        steering;
        label = "Jerusalem cross Wang-Werner",
        options,
        seed,
        run_direct,
    )
end


results = run(run_direct = haskey(ENV, "RUN_DIRECT"));
#results = run();
plot_comparison(results; i = 4, j = 2, savepath = joinpath(@__DIR__, "jerusalem_wang_werner.png"))
