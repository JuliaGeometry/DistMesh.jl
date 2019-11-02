
function plotout(statsdata, qualities, name::AbstactString)
    timestamp = Dates.format(now(), "yyyy-mm-ddTHH_MM")

    qual_hist = Plots.histogram(qualities, title = "Quality", bins=30, legend=false)
    avg_plt = Plots.plot(statsdata.average_qual, title = "Average Quality", legend=false, ylabel="Quality")
    vline!(statsdata.retriangulations, line=(0.2, :dot, [:red]))
    med_plt = Plots.plot(statsdata.median_qual, title = "Median Quality", legend=false, ylabel="Quality")
    vline!(statsdata.retriangulations, line=(0.2, :dot, [:red]))
    maxdp_plt = Plots.plot(statsdata.maxdp, title = "Max Displacement", legend=false, ylabel="Edge Displacement")
    vline!(statsdata.retriangulations, line=(0.2, :dot, [:red]))
    maxmove_plt = Plots.plot(statsdata.maxmove, title = "Max Move", legend=false, ylabel="Point Displacement")
    vline!(statsdata.retriangulations, line=(0.2, :dot, [:red]))
    plt = Plots.plot(avg_plt, med_plt,maxdp_plt,maxmove_plt,layout=(2,2), xlabel="Iteration")

    # TODO localize to dir
    savefig(plt, "output/result_stat$name$timestamp.svg")
    savefig(qual_hist, "output/result_qualmaxmove$name$timestamp.svg")
end