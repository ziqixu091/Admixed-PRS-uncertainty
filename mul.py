import pandas as pd
import xarray as xr
import numpy as np

# 528612 snp, 500 individuals
raw_beta = pd.read_csv('beta.tsv.gz', sep='\t')

raw_geno = xr.open_zarr("ukb_eur_afr.zip")
#raw_geno = admix.data.load_lab_dataset("ukb_eur_afr")

#raw_geno = dset.to_dataframe()

snp_info = raw_geno["snp"].values  # 'chr1:831716:A:G'
r_hap = raw_geno["geno"].values # 4327 ind, 528632 snp, 2 hap
r_n_ind = np.shape(r_hap)[0]
r_n_snp = np.shape(r_hap)[1]
r_geno = r_hap[:, :, 0] + r_hap[:, :, 1] # 4327 ind, 528632 snp
r_geno = r_geno.transpose() # snp: 528632, indiv: 4327
r_geno = r_geno.values
# r_geno = np.zeros((r_n_ind, r_n_snp)) 
# for x in range(0, r_n_ind-1):
#     for y in range(0, r_n_snp-1):
#         r_geno[x][y] = hap[x][y][0] + hap[x][y][1]


# raw beta and raw geno needs to be refined
# before matching, they should be in form of:
#    geno (n_ind, n_snp)  (500 out of 4327, 528612 out of 528632)  
#    beta (n_ind, n_snp)  (500, 528612)

# initialization 
n_indiv = 500
beta = [0]*500 # (snp, 500)
geno = [0]*500 # (snp, 500)

# do refinement of SNPs by checking if the alleles matches
for i in range(1, 22):
    for raw_beta.CHR == i:
        for item in raw_beta.SNP:
            b_r_i = raw_beta.index[(raw_beta.SNP==item)==True]
            a1 = raw_beta["A1"][b_r_i].to_string()
            a2 = raw_beta["A2"][b_r_i].to_string()
            alleles = [a1[5], a2[5]]
            
            line = snp_info[b_r_i[0]]
            if line[3] == i:
                if line[-3] in alleles and line[-1] in alleles: 
                    row = raw_beta.loc[raw_beta["SNP"] == item].values.tolist()[0]
                    # beta = np.append(beta, row[5:], axis=0)
                    beta = np.vstack([beta, row[5:]])
                    
                    row_2 = r_geno[b_r_i[0]]
                    geno = np.vstack([geno, row_2[:500]])
                    

beta = np.delete(beta, 0, axis=0)
geno = np.delete(geno, 0, axis=0)

# after refinement, beta and geno should be in the same shape

heritability = 0.5
phe_g = beta*geno
phe_e = np.zeros_like(phe_g) 

var_g = np.var(phe_g[:, 1])
var_e = var_g * ((1.0 / heritability) - 1)
phe_e[:, 1] = np.random.normal(loc=0.0, scale=np.sqrt(var_e), size=n_indiv) 

pheno = phe_g + phe_e


# estimation with real phenotype

