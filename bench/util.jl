
function plotout(statsdata, qualities, folder, name)

    qual_hist = Plots.histogram(qualities, title = "Quality", bins=30, legend=false)
    avg_plt = Plots.plot(statsdata.average_qual, title = "Average Tri Quality", legend=false, ylabel="Quality")
    vline!(statsdata.retriangulations, line=(0.2, :dot, [:red]))

    med_plt = Plots.plot(statsdata.median_qual, title = "Median Tri Quality", legend=false, ylabel="Quality")
    vline!(statsdata.retriangulations, line=(0.2, :dot, [:red]))

    min_plt = Plots.plot(statsdata.minimum_qual, title = "Minimum Tri Quality", legend=false, ylabel="Quality")
    vline!(statsdata.retriangulations, line=(0.2, :dot, [:red]))

    max_plt = Plots.plot(statsdata.maximum_qual, title = "Maximum Tri Quality", legend=false, ylabel="Quality")
    vline!(statsdata.retriangulations, line=(0.2, :dot, [:red]))

    min_tq_plt = Plots.plot(statsdata.min_volume_edge_ratio, title = "Minimum Tet Quality", legend=false, ylabel="Vol/Edge Ratio")
    vline!(statsdata.retriangulations, line=(0.2, :dot, [:red]))

    max_tq_plt = Plots.plot(statsdata.max_volume_edge_ratio, title = "Maximum Tet Quality", legend=false, ylabel="Vol/Edge Ratio")
    vline!(statsdata.retriangulations, line=(0.2, :dot, [:red]))

    maxdp_plt = Plots.plot(statsdata.maxdp, title = "Max Displacement", legend=false, ylabel="Edge Displacement")
    vline!(statsdata.retriangulations, line=(0.2, :dot, [:red]))

    maxmove_plt = Plots.plot(statsdata.maxmove, title = "Max Move", legend=false, ylabel="Point Displacement")
    vline!(statsdata.retriangulations, line=(0.2, :dot, [:red]))

    plt = Plots.plot(avg_plt, med_plt,min_plt,max_plt,min_tq_plt,max_tq_plt,maxdp_plt,maxmove_plt,layout=(4,2), xlabel="Iteration")

    savefig(plt, "$folder/result_stat$name.svg")
    savefig(qual_hist, "$folder/result_qual$name.svg")
end
