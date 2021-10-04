from os.path import join
import xarray as xr
import numpy as np
import pandas as pd
import admix
import fire
np.random.seed(1234)


# read in command line and simulate
# python3 simulate.py --data_dir ${PLINK_DIR} --chr_i ${chr_i} --var_g ${var_g} --var_e ${var_e}
def simulate(data_dir, chr_i, var_g, var_e):
    '''
    print("read in data dir: ", data_dir)
    print("read in chr num: ", chr_i)
    print("read in var g: ", var_g)
    print("read in var e: ", var_e)
    '''
    # read in plink data using pandas_plink and set up xarray Dataset
    dset_eur_train = admix.io.read_plink(data_dir+"/eur_train/chr"+str(chr_i))
    dset_eur_val = admix.io.read_plink(data_dir+"/eur_val/chr"+str(chr_i))
    dset_eur_test = admix.io.read_plink(data_dir+"/eur_test/chr"+str(chr_i))
    dset_admix = admix.io.read_plink(data_dir+"/admix/chr"+str(chr_i))

    # add up 
    dset = xr.concat([dset_eur_train, dset_eur_val, dset_eur_test, dset_admix], dim="indiv")

    # simulate phenotype
    sim = admix.simulate.continuous_pheno_1pop(
        dset, var_g=var_g, var_e=var_e, n_sim=10
    )

    beta = sim["beta"].copy()
    pheno = sim["pheno"].copy()
    pheno.index = pheno.index.str.split("_", expand=True)
    pheno.index.names = ["FID", "IID"]
    pheno = pheno.reset_index()

    beta.to_csv("beta.csv")
    pheno.to_csv("pheno.csv")
    
    admix.tools.plink_assoc(
        data_dir+"/eur_train/chr"+str(chr_i), pheno, out_prefix="out/plink"
    )


if __name__ == '__main__':
  fire.Fire(simulate)


