import pandas as pd
import gseapy as gs


def process_output(counts, signature, phenodata2, phenodata1, outfilename): 
    """
    counts:: file of for expression
    signature:: signature file in mt format
    pheno:: pheno specific file
    outpfile:: what to save outfile at
    """
    gene_list = pd.read_csv(counts, index_col=0).T
    sig = pd.read_csv(signature, skiprows=1)
    sig.columns = ["sig"]
    sig_dict = {"cluster6": sig["sig"].to_list()}
    ss = gs.ssgsea(data=gene_list,
                   gene_sets=sig_dict,
                   outdir=None,
                   sample_norm_method='custom', # choose 'custom' will only use the raw value of `data`
                   no_plot=True)

    pheno = pd.read_csv(phenodata1, index_col=0)
    pheno["response"] = pheno["pd_6m"].map({"yes":"NR", "no":"R"})
    pheno["enrichment_C6"] = ss.res2d.T["cluster6"]
    pheno["bin_surv"] = pheno["enrichment_C6"] > pheno["enrichment_C6"].median()


    survival_data = pd.read_csv(phenodata2, index_col=0)
    survival_data["patient.1"] = survival_data["patient.2"]
    merged_df = pheno.merge(survival_data)
    condensed = merged_df[["patient.1", "pd_6m", "response", "enrichment_C6", "bin_surv", "adjir_status", "adj_time_ir", "os_status", "os_time"]]
    condensed.to_csv(outfilename)
    

process_output("data/timepoint0.csv", 
               "data/pace.grp", 
               "data/phenodata_all.csv", 
               "data/phenodata_PLS.csv", "results/pace_Timepoint1.csv")

process_output("data/timepoint1.csv", 
               "data/pace.grp", 
               "data/phenodata_all.csv", 
               "data/phenodata_PLS.csv", "results/pace_Timepoint0.csv")

process_output("data/timepoint2.csv", 
               "data/pace.grp", 
               "data/phenodata_all.csv", 
               "data/phenodata_PLS.csv", "results/pace_Timepoint2.csv")