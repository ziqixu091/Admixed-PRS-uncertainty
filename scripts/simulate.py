import os
import xarray as xr
import numpy as np
import pandas as pd
from tqdm import tqdm
from os.path import join
import admix
import fire
np.random.seed(1234)

# helper function to transfrom PLINK to ldpred2
def plink2_assoc_to_ldpred2(plink2_assoc_path):
    assoc = pd.read_csv(plink2_assoc_path, delim_whitespace=True)
    assert np.all(assoc["A1"] == assoc["ALT"])
    assoc = assoc[
        [
            "#CHROM",
            "ID",
            "POS",
            "ALT",
            "REF",
            "OBS_CT",
            "BETA",
            "SE",
            "T_STAT",
            "P",
        ]
    ].rename(
        columns={
            "#CHROM": "CHR",
            "ID": "SNP",
            "POS": "BP",
            "ALT": "A1",
            "REF": "A2",
            "OBS_CT": "N",
            "T_STAT": "STAT",
        }
    )
    return assoc


# read in command line and simulate
# python3 simulate.py --data_dir ${PLINK_DIR} --chr_i ${chr_i} --h2 ${} --cau_prop ${} --out {$dir}
# 10% h2 fix var_e to 1, 0.1, 0 to 1, h2 = var_g/(var_g+var_e)
# cau_prop = cau / n_snp, 0.01
def simulate(data_dir, chr_i, h2, cau_prop, out):
    
    print("read in data dir: ", data_dir)
    print("read in chr num: ", chr_i)
    print("read in h2: ", h2)
    print("read in cau_prop: ", cau_prop)
    
    # read in plink data using pandas_plink and set up xarray Dataset
    dset_eur_train = admix.io.read_plink(data_dir+"/eur_train/chr"+str(chr_i))
    dset_eur_val = admix.io.read_plink(data_dir+"/eur_val/chr"+str(chr_i))
    dset_eur_test = admix.io.read_plink(data_dir+"/eur_test/chr"+str(chr_i))
    dset_admix = admix.io.read_plink(data_dir+"/admix/chr"+str(chr_i))
    
    # add up 
    dset = xr.concat([dset_eur_train, dset_eur_val, dset_eur_test, dset_admix], dim="indiv")

    # derive var_g, var_e, n_causal
    n_snp = dset_eur_train["snp"].shape[0]
    n_causal = round(n_snp * cau_prop)
    
    # fix var_e to 1
    var_e = 1
    var_g = h2/(1-h2)
    
    # simulate phenotype
    sim = admix.simulate.continuous_pheno_1pop(
        dset, var_g=var_g, var_e=var_e, n_sim=10, n_causal = n_causal
    )
    
    # output simulation results
    beta = sim["beta"].copy()
    pheno = sim["pheno"].copy()
    pheno.index = pheno.index.str.split("_", expand=True)
    pheno.index.names = ["FID", "IID"]
    pheno = pheno.reset_index()
    
    outdir = str(out)
    beta_out = 'beta_h2_'+str(h2)+'.csv'
    pheno_out = 'pheno_h2_'+str(h2)+'.csv'
    plink_out = outdir+'plink'
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    str_beta = os.path.join(outdir, beta_out)    
    str_pheno = os.path.join(outdir, pheno_out)
    
    beta.to_csv(str_beta)
    pheno.to_csv(str_pheno)
    
    admix.tools.plink_assoc(
        data_dir+"/eur_train/chr"+str(chr_i), pheno, out_prefix=plink_out
    )
    
    # formating the plink data for ldpred2
    for i_sim in range(10):
        path = str_plink+'SIM_'+str(i_sim)+'.glm.linear'
        assoc_ldpred2 = plink2_assoc_to_ldpred2(path)
        str_assoc = outdir+'assoc.sim'+str(i_sim)+'.ldpred2.csv'
        assoc_ldpred2.to_csv(str_assoc)
    

if __name__ == '__main__':
  fire.Fire(simulate)



