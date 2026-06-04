include(joinpath(@__DIR__, "common.jl"))
using PSSFSS.Constants: μ₀

function cpss_meander_strip_case(; nfreq = 101, ntri = 600)
    P = 5.2
    d1 = 2.61
    d2 = 3.81
    h0 = 2.44
    h2 = 2.83
    w0x = 0.46
    w0y = 0.58
    w1 = 0.21
    w2x = 0.25
    w2y = 0.17
    σ = 58e6
    units = mm

    outer(orient) = meander(;
        a = P,
        b = P,
        w1 = w2y,
        w2 = w2x,
        h = h2 + w2x,
        units,
        ntri,
        σ,
        orient,
    )
    inner = meander(;
        a = P,
        b = P,
        w1 = w0y,
        w2 = w0x,
        h = h0 + w0x,
        units,
        ntri,
        σ,
    )
    strip(orient) = diagstrip(; P, w = w1, units, Nl = 60, Nw = 4, orient, σ)

    substrate = Layer(width = 0.127mm, epsr = 2.17, tandel = 0.0009)
    foam(w) = Layer(width = w, epsr = 1.043, tandel = 0.0017)
    sheets = [outer(-90), strip(-45), inner, strip(45), outer(90)]

    strata = [
        Layer()
        substrate
        sheets[1]
        foam(d2 * 1mm)
        substrate
        sheets[2]
        foam(d1 * 1mm)
        sheets[3]
        substrate
        foam(d1 * 1mm)
        substrate
        sheets[4]
        foam(d2 * 1mm)
        sheets[5]
        substrate
        Layer()
    ]

    return strata, collect(range(10.0, 20.0, length = nfreq)), (θ = 0, ϕ = 0)
end

function run(; nfreq = 101, ntri = 600, seed = 1, run_direct = false, options = DEFAULT_SWEEP_OPTIONS)
    strata, freqs, steering = cpss_meander_strip_case(; nfreq, ntri)
    return run_sweep_comparison(
        strata,
        freqs,
        steering;
        label = "Meander/strip CPSS",
        options,
        seed,
        run_direct,
    )
end

results = run();
plot_comparison(results; savepath = joinpath(@__DIR__, "cpss_meander_strip.png"))
