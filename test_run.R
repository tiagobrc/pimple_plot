source("example_gen.R")
source("pimple_plot.R")

plot_volcano(df = as.data.frame(resLFC), padj.sig = 0.05, avg_fc.sig = 0.5,
             name = "condition_treated_vs_untreated",
             highlight_genes = NULL, highlight_top = 10,
             seed = 16)


