# RA_joint
Leveraging Hardy-Weinberg disequilibrium to test association in case-control studies

## Table of Contents
- [Import RA joint source functions](#import_functions)
- [Testing association with RA joint](#RA_joint)

## :heavy_plus_sign: Examples
- <a name="import_functions"></a> Save [RA joint functions](https://github.com/lzhangdc/RA_joint/blob/main/RA_joint.R) in the working directory
```R
source('RA_joint.R')
```

- <a name="indep"></a> [Independent samples], Y can be continuous or binary
```R
########################################################################
########################################################################
#### Function inputs:
#### - g: a vector of genotypes taking values 0, 1, 2 
#### - y: a binary vector of phenotypes
########################################################################
########################################################################
#### Function output:
#### - RA_joint: the p-value by proposed RA joint test 
#### - genotypic: genotypic Wald's test p-value 
#### - geno_score: genotypic score test p-value 
#### - add_wald: additive Wald's test p-value
#### - add_score: additive score test p-value
#### - control_maf: maf estimated in the control sample
#### - control_delta: HWD measure in the control sample
#### - control_HWE: p-value of HWE test in the control sample only
#### - case_maf: maf estimated in the case sample
#### - case_delta: HWD measure in the case sample
#### - case_HWE: p-value of HWE test in the case sample only
#### - pool_maf: maf estimated in the pooled sample
#### - pool_delta: HWD measure in the pooled sample
#### - pool_HWE: p-value of HWE test in the pooled sample
########################################################################

temp_result <- fast_MI_joint_HWE(gSNP=g, pheno=y)
names(temp_result) <- c('RA_joint', 'genotypic', 'geno_score', 'add_wald', 'add_score', 'control_maf','control_delta', 'control_HWE', 'case_maf','case_delta', 'case_HWE','pool_maf','pool_delta', 'pool_HWE')

```
