using DelimitedFiles: writedlm
function tcase_not_pssfss(str::AbstractString)
    strs = split(str)
    for (i, s) in pairs(strs)
        strs[i] = s == "PSSFSS" ? s : titlecase(s)
    end
    return join(strs, " ")
end

function make_summarycsv(fname::AbstractString)
    (dir, _) = splitdir(fname)
    outfile = joinpath(dir, "summary.csv")

    # Determine examples that were run
    examples = String[]
    for l in eachline(fname)
        occursin("Example:", l) || continue
        push!(examples, split(l, ":") |> last |> strip)
    end

    headings = ["Frequency grid",
                "Loewner total time", "Random total time", "Direct grid time", "PSSFSS fast time",
                "Loewner samples", "Random samples", "PSSFSS fast samples",
                "Loewner max relerr", "Loewner med relerr",
                "Random max relerr", "Random med relerr",
                "PSSFSS fast maxerr", "PSSFSS fast mederr",
                ]

    results = Matrix{Any}(undef, 1 + length(examples), 1 + length(headings))
    results[1,1] = "Example Case"
    results[1, 2:end] .= tcase_not_pssfss.(headings)
    results[2:end, 1] .= tcase_not_pssfss.(examples)


    open(fname, "r") do fid
        line = ""
        for (rowm1, ex) in pairs(examples)
            while !occursin(ex, line); (line = readline(fid)); end
            for (colm1, heading) in pairs(headings)
                while !occursin(heading, line)
                    line = readline(fid)
                end
                xstr = split(line, ":")[2] |> split |> first
                T = occursin(".", xstr) ? Float64 : Int
                results[1 + rowm1, 1 + colm1] = parse(T, xstr)
            end
        end
    end

    writedlm(outfile, results, ",")
    return results
end