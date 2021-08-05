import pandas as pd
import numpy as np
import xarray as xr

# manage space 
from dask.distributed import Client
import dask
client = Client(processes=False, threads_per_worker=1, n_workers=1, memory_limit='16GB')


# read in 528612 snp, 500 samples
raw_beta = pd.read_csv("beta.tsv.gz", sep="\t")
raw_geno = xr.open_zarr("ukb_eur_afr.zip")

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

# remove nan values in real height data 
height = geno_aligned["height@indiv"].values # 4327 indiv
nan_li = np.argwhere(~np.isnan(height))
height = height[~np.isnan(height)] # 4296 indiv
nan_li = np.concatenate(nan_li)
geno = geno[nan_li, :, :] # (4296, 511278, 2)

# phenotype: [4296 indiv, 500 samples]
# correlation: [4296 indiv, 500 samples]
pheno = np.matmul(geno.sum(axis=2), beta_aligned[f"SAMPLE_{i}" for i in range (1, 501)].values).compute()  
corr = np.corrcoef(height, pheno[:, sample_i]) for sample_i in range(500)

corr.to_csv("corr.csv")


# plot

# more traits 





