import pandas as pd
import numpy as np
import xarray as xr
from dask.distributed import Client
import dask
import fire


# read in command line and process
def process(trait, beta, geno, out):
    """
    print("read in beta file: ", beta)
    print("read in geno file: ", geno)
    print("read in output file: ", out)
    print("reaad in trait: ", trait)
    """

    # manage space 
    client = Client(processes=False, threads_per_worker=1, n_workers=1, memory_limit='28GB')


    # read in 528612 snp, 500 samples of geno, beta
    raw_beta = pd.read_csv(beta, sep="\t")
    raw_geno = xr.open_zarr(geno)

    # transofrm into dataframe with corresponding columns
    geno_info = raw_geno.snp.to_dataframe().rename(
        columns={"a0@snp": "A1", "a1@snp": "A2", "chrom@snp": "CHR", "pos_hg19@snp": "POS"}
    )
    # extract the snp info
    beta_info = raw_beta[["CHR", "SNP", "POS", "A1", "A2"]]
    # merge do the match job
    snp_merged = pd.merge(beta_info, geno_info, on=["CHR", "POS", "A1", "A2"])

    # refinement
    id_beta = snp_merged["SNP"].values # rs12562034
    id_geno = snp_merged["snp"].values # chr1:833068:G:A

    # generate refined beta and geno data
    beta_aligned = raw_beta[raw_beta.SNP.isin(id_beta)]  # 511278 snp, 500 samples
    geno_aligned = raw_geno.sel(snp=id_geno)  # indiv: 4327 ploidy: 2 snp: 511278

    # casting data type 
    geno = geno_aligned.geno.data.astype(float)
    lanc = geno_aligned.lanc.data

    # remove nan values indiv in real height data and lanc
    trait = geno_aligned[f"{trait}@indiv"].values # 4327 indiv
    nan_li = np.argwhere(~np.isnan(trait))
    trait = trait[~np.isnan(trait)] # 4296 indiv
    nan_li = np.concatenate(nan_li)
    geno = geno[nan_li, :, :] # (4296, 511278, 2) warning? PerformanceWarning: Slicing is producing a large chunk. To accept the large chunk and silence this warning, set the option
    lanc = lanc[nan_li, :, :] # (4296, 511278, 2)

    # phenotype: [4296 indiv, 500 samples]
    # correlation: [4296 indiv, 500 samples]
    pheno = np.matmul(geno.sum(axis=2), beta_aligned[f"SAMPLE_{i}" for i in range (1, 501)].values).compute()  # bug wrong for
    pheno = []
    for i in range(1, 500):
        one = np.matmul(geno.sum(axis=2), beta_aligned[f"SAMPLE_{i}"].values).compute()
        pheno.append(one)
    corr = np.corrcoef(height, pheno[:, sample_i]) for sample_i in range(1, 501) 

    np.mean(corr)
    np.std(corr)
    # For example:
    # height: correlation: 0.11, std: 0.01
    # LDL: correlation: 0.11, std: 0.01
    # HDL:
    # Or a figure of histogram:
    plt.hist(corr)

    # corr.to_csv("corr.csv")

if __name__ == '__main__':
    fire.Fire(process)

# python new.py --trait height --beta beta.tsv.gz --geno ukb_eur_afr.zip --out trait_corr.csv


# plot
# matplotlib / seaborn
# Things to plot
# 1. Scatter plot:
# x-axis: height (true height in dataset),
# y-axis: prediction +- 2 * standard deviation

# correlation follows a normal distribution?? +-2sd  95%?
# calculate standard deviation 

import matplotlib.pyplot as plt
plt.scatter(height, pheno, facecolor="none", edgecolor="b", s=50, label="first plot")
plt.plot(height, pheno.mean, c="r", label="mean")
plt.fill_between(height, pheno - 2*pheno.std, pheno + 2*pheno.std, color="pink", label="2std.", alpha=0.5)
plt.legend()
plt.show()

# 2. Scatter plot
# x-axis: average local ancestry (the proportion of african ancestry in each individual) 0/1
# y-axis: height
# Question: 
hap = lanc.sum(axis=2)
avg_lanc = hap.sum(axis=1)/(511278*2)
plt.scatter(avg_lanc, height, facecolor="none", edgecolor="b", s=50, label="second plot")
plt.plot(avg_lanc, height.mean, c="r", label="mean")
plt.fill_between(avg_lanc, height - 2*height.std, height + 2*height.std, color="pink", label="2std.", alpha=0.5)
plt.legend()
plt.show()

# 3. Scatter plot
# x-axis: average local ancestry (the proportion of african ancestry in each individual)
# y-axis: predicted height +- 2 * standard deviation

plt.scatter(avg_lanc, pheno, facecolor="none", edgecolor="b", s=50, label="third plot")
plt.plot(avg_lanc, pheno.mean, c="r", label="mean")
plt.fill_between(avg_lanc, pheno - 2*pheno.std, pheno + 2*pheno.std, color="pink", label="2std.", alpha=0.5)
plt.legend()
plt.show()

# why no corr?








