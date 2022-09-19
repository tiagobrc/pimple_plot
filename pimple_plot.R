plot_volcano <- function(df, name, highlight_genes = NULL, highlight_num = 10, highlight_top = T, seed = 42, padj.sig = 0.1, avg_fc.sig = 0){
  
  require("ggplot2")
  require("dplyr")
  require("ggrepel")
  
  ggplot.theme <- theme(axis.text=element_text(size=16, family = "Helvetica", colour = "black"),
                        axis.title=element_text(size=16, family = "Helvetica", colour = "black"),
                        legend.title = element_text(size=16, colour = "black"),
                        legend.text=element_text(size=16, colour = "black"),
                        plot.title = element_text(size=16, family = "Helvetica"),
                        panel.grid = element_line(colour = "grey92"),
                        panel.background = element_rect(fill = "white", colour = NA))
  
  df.1 <- df %>% filter(log2FoldChange > avg_fc.sig & padj < padj.sig)
  df.2 <- df %>% filter(log2FoldChange < avg_fc.sig*-1 & padj < padj.sig)
  name1 <- gsub(paste0(gsub("_.*","",name), "_(.*)_vs.*"), "\\1", name)
  name2 <- gsub(paste0(gsub("_.*","",name), "_.*_vs_(.*)"), "\\1", name)
  name.y <- max(-log10(df$padj), na.rm = T) * 0.8
  name.x1 <- range(df$log2FoldChange)[2] * 0.5
  name.x2 <- range(df$log2FoldChange)[1] * 0.5
  
  if(!is.null(highlight_genes)){
    genes <- highlight_genes
    df$genes <- rownames(df)
    df.sub <- df[genes,]
  }else{
    if(highlight_top){
      df$genes <- rownames(df)
      df.sub <- df %>% filter(abs(log2FoldChange) > 0 & padj < 0.1) 
      
      df.sub <- df.sub %>% slice_min(n = highlight_num, order_by = padj)
      
    }else{
      df$genes <- rownames(df)
      df.sub <- df %>% filter(abs(log2FoldChange) > 0 & padj < 0.1)
      df.sub <- df.sub %>% arrange(padj)
      set.seed(seed)
      df.sub <- df.sub[sample(1:nrow(df.sub), sample(1:highlight_num, 1, replace = F), replace = F),]
      }
    }
  
  p <- ggplot(data = df, name = name) +
    geom_point(df, mapping = aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(df.1, mapping = aes(x = log2FoldChange, y = -log10(padj)), color = "red2") +
    geom_point(df.2, mapping = aes(x = log2FoldChange, y = -log10(padj)), color = "midnightblue") +
    annotate("text", x = name.x1, y = name.y, label = c(paste0(nrow(df.1), " ASVs high in ", name1)), size = 4, color = "red2") + 
    annotate("text", x = name.x2, y = name.y, label = c(paste0(nrow(df.2), " ASVs high in ", name2)), size = 4, color = "midnightblue")  +
    geom_text_repel(df.sub, mapping = aes(x = log2FoldChange, y = -log10(padj)), label = df.sub$genes, 
                    point.padding = 0.2,
                    nudge_x = .15,
                    nudge_y = .5,
                    segment.curvature = -1e-20,
                    arrow = arrow(length = unit(0.015, "npc")), box.padding = 1, max.overlaps = Inf) +
    ggplot.theme +
    theme_bw()
  
  return(p)
}