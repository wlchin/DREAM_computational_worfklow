
library(ComplexHeatmap)

C2_C3_r = read.csv("DE_results/C2_C3_r.csv")
Base_C2_r = read.csv("DE_results/Base_C2_r.csv")
Base_C2_nr = read.csv("DE_results/Base_C2_nr.csv")
C2_C3_nr = read.csv("DE_results/C3_C2_nr.csv")
Base_r_nr = read.csv("DE_results/Base_r_vs_nr.csv")
C2_r_vs_nr = read.csv("DE_results/C2_C3_r.csv")
C3_r_vs_nr = read.csv("DE_results/C3_r_vs_nr.csv")

filtered_output <- function(x){
  filtered_df <- x[(abs(x$logFC) > 1.0) & (x$adj.P.Val < 0.05),]$X
  filtered_df
}

lt = list(R.Timepoint1.vs.Timepoint2 = filtered_output(C2_C3_r),
          R.Timepoint0.vs.Timepoint1 = filtered_output(Base_C2_nr),
          NR.Timepoint0.vs.Timepoint1 = filtered_output(Base_C2_r),
          NR.Timepoint1.vs.Timepoint2 = filtered_output(C2_C3_nr),
          R.vs.NR.Timepoint1 = filtered_output(C2_r_vs_nr),
          R.vs.NR.Timepoint2 = filtered_output(C3_r_vs_nr),
          R.vs.NR.Timepoint0 = filtered_output(Base_r_nr))
          
list_to_matrix(lt)

m1 = make_comb_mat(lt)

png("results/upsetplot.png", width = 4.8, height = 4.27, units = "in", res = 600)
UpSet(m1)
dev.off()

