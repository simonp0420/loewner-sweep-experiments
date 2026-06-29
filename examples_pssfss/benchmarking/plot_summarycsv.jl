using DelimitedFiles: readdlm
using GLMakie

function get_settings(fname)
    sweep_settings = ""
    julia_version = ""
    la_lib = ""
    j_nthreads = ""
    cpu = ""
    for line in eachline(fname)
        sline = strip(line)
        isempty(sline) && continue

        if sline[1] == '('
            sweep_settings = sline
            continue
        elseif occursin("Julia Version", sline)
            julia_version = sline
            continue
        elseif occursin("BLAS", sline)
            la_lib = "BLAS: Unknown"
            occursin("mkl", sline) && (la_lib = "BLAS: MKL")
            occursin("openblas", sline) && (la_lib = "BLAS: OpenBLAS")
            continue
        elseif occursin("Threads:", sline)
            j_nthreads = "Threads.nthreads = " * split(sline)[2]
            continue
        elseif occursin("CPU: ", sline)
            cpu = sline
            continue
        end

        !isempty(sweep_settings) && !isempty(julia_version) && !isempty(la_lib) &&
        !isempty(j_nthreads) && !isempty(cpu) && break
    end

    label = join((cpu, julia_version, j_nthreads, la_lib), ", ") * "\n" * sweep_settings
    return label
end

function plot_summarycsv(fname::AbstractString)
    fname = abspath(fname)
    (dir, _) = splitdir(fname)
    lastdir = last(splitpath(dir))
    plotfile = joinpath(dir, lastdir * ".png")
    summaryfile = joinpath(dir, "runall_results.txt")
    settings = get_settings(summaryfile)
    srch_str = "benchmarking"
    i = findfirst(srch_str, plotfile)[end] + 2
    settings = string(replace(plotfile[i:end], "\\" => "/"), '\n', settings)

    mat = readdlm(fname, ',')

    fig = Figure(size = (1000, 620))

    # Errors wrt PSSFSS Discrete Sweep
    loewnermaxerr_col = findfirst(==("Loewner Max Relerr"), @view mat[1,:])
    randommaxerr_col = findfirst(==("Random Max Relerr"), @view mat[1,:])
    pssfssmaxerr_col = findfirst(==("PSSFSS Fast Maxerr"), @view mat[1,:])

    loewner_maxerr = @views mat[2:end, loewnermaxerr_col]
    random_maxerr = @views mat[2:end, randommaxerr_col]
    pssfss_maxerr = @views mat[2:end, pssfssmaxerr_col]

    ax1 = Axis(
        fig[1, 1];
        xticks = (1:length(loewner_maxerr), ["" for _ in 1:length(loewner_maxerr)]),
        ylabel = "log₁₀(error)",
        yticks = -10:2:10,
        title = """Max Relerr for Case "$(lastdir)\"""",
    )

    scatter!(ax1, 1:length(loewner_maxerr), log10.(loewner_maxerr);
        marker = :rect,
        markersize = 14,
        color = :red,
        strokewidth = 0.8,
        strokecolor = :white,
        label = "Loewner")

    scatter!(ax1, 1:length(random_maxerr), log10.(random_maxerr);
        marker = :circle,
        markersize = 12,
        color = :blue,
        strokewidth = 0.8,
        strokecolor = :white,
        label = "Random")

    scatter!(ax1, 1:length(pssfss_maxerr), log10.(pssfss_maxerr);
        marker = :utriangle,
        markersize = 10,
        color = :green,
        strokewidth = 0.8,
        strokecolor = :white,
        label = "PSSFSS Fast")

    #axislegend(ax1, position = :lt)
    Legend(fig[1,2], ax1, "", halign = :left)

    # Timing relative to PSSFSS Fast Sweep
    loewnertime_col = findfirst(==("Loewner Total Time"), @view mat[1,:])
    randomtime_col = findfirst(==("Random Total Time"), @view mat[1,:])
    pssfsstime_col = findfirst(==("PSSFSS Fast Time"), @view mat[1,:])

    loewner_rel_time = @views mat[2:end, loewnertime_col] ./ mat[2:end, pssfsstime_col]
    random_rel_time = @views mat[2:end, randomtime_col] ./ mat[2:end, pssfsstime_col]

    ax2 = Axis(
        fig[2, 1];
        ylabel = "Normalized T_tot",
        title = """Timing for Case "$(lastdir)" Relative to PSSFSS Fast Sweep""",
        xticks = (1:length(loewner_rel_time), mat[2:end, 1]),
        xticklabelrotation = -3π/12,
    )

    scatter!(ax2, 1:length(loewner_rel_time), loewner_rel_time;
        marker = :rect,
        markersize = 14,
        color = :red,
        strokewidth = 0.8,
        strokecolor = :white,
        label = "Loewner")

    scatter!(ax2, 1:length(random_rel_time), random_rel_time;
        marker = :circle,
        markersize = 12,
        color = :blue,
        strokewidth = 0.8,
        strokecolor = :white,
        label = "Random")

    hlines!(ax2, 1.0, color = :gray, linestyle = :dash, linewidth = 4)
    #axislegend(ax2, position = :rt)
    Legend(fig[2,2], ax2, "", halign = :left)


    # Settings:
    Label(fig[3,:], settings, fontsize = 12, halign = :center)
    colgap!(fig.layout, 12)

    display(fig)
    save(plotfile, fig)
    println(plotfile)

    return(fig)
end

function plot_all_summarycsvs(startdir = pwd())
    searchname = "summary.csv"
    for (root, dirs, files) in walkdir(startdir)
        if searchname in files
            summary_file = joinpath(root, searchname)
            plot_summarycsv(summary_file)
        end
    end
end