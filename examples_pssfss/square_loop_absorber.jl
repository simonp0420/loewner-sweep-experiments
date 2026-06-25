include(joinpath(@__DIR__, "common.jl"))

function square_loop_absorber_case(; design = :thin, nfreq = 121, ntri = 750)
    designs = Dict(
        :thin => (i = 1, Rs = 15.0),
        :medium => (i = 2, Rs = 40.0),
        :thick => (i = 3, Rs = 70.0),
    )
    haskey(designs, design) || error("design must be one of :thin, :medium, :thick")

    D = 11.0
    idx = designs[design].i
    r_outer = sqrt(2) / 2 * D / 8 * [5, 6, 7]
    thickness = D / 16 * [1, 2, 3]
    r_inner = r_outer - sqrt(2) .* thickness

    sheet = polyring(;
        sides = 4,
        s1 = [D, 0],
        s2 = [0, D],
        ntri,
        orient = 45,
        a = [r_inner[idx]],
        b = [r_outer[idx]],
        Zsheet = designs[design].Rs,
        units = mm,
    )

    strata = [
        Layer()
        sheet
        Layer(width = 5mm)
        pecsheet()
        Layer()
    ]

    return strata, collect(range(1.0, 25.0, length = nfreq)), (ϕ = 0, θ = 0)
end

function run(; design = :thin, nfreq = 121, ntri = 750, seed = 1, run_direct = false, options = DEFAULT_SWEEP_OPTIONS)
    strata, freqs, steering = square_loop_absorber_case(; design, nfreq, ntri)
    return run_sweep_comparison(
        strata,
        freqs,
        steering;
        label = "Square loop absorber $(design)",
        options,
        seed,
        run_direct,
    )
end

results = run(run_direct = haskey(ENV, "RUN_DIRECT"));
#results = run();
plot_comparison(results; savepath = joinpath(@__DIR__, "square_loop_absorber.png"))
