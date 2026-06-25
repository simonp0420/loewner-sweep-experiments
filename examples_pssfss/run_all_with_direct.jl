cd(@__DIR__)
using Pkg: Pkg
Pkg.activate(pwd())

using ShareAdd: @usingany
@usingany MKL

println(versioninfo())
using LinearAlgebra: BLAS
println("  BLAS: $(BLAS.get_config())\n")

to_namedtuple(obj) = NamedTuple{propertynames(obj)}(getproperty.(Ref(obj), propertynames(obj)))
include(joinpath(@__DIR__, "common.jl")) # To define DEFAULT_SWEEP_OPTIONS

to_namedtuple(DEFAULT_SWEEP_OPTIONS) |> println
println()
withenv("RUN_DIRECT" => true) do
    include("band_pass_filter.jl")
    include("cpss_meander_strip.jl")
    include("cpss_meander.jl")
    include("cross_substrate.jl")
    include("flexible_absorber.jl")
    include("jerusalem_wang_werner.jl")
    include("resistive_square_patch.jl")
    include("splitring.jl")
    include("square_loop_absorber.jl")
    include("symmetric_strip.jl")
end