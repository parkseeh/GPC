## GPC
An R package for calculating the GWAS power based on the ddds ratio, minor allele frequency and number of sample.

## Intallation
To install directly from Github:
```r
require(devtools)
devtools::install_github("CESP-ExpHer/GPC")
```

## Description
Please read the description using 'help()' in R 

```r
library(GPC)
help(GPC)
```


## Example
#### 1. Linear Model
```r
# If the GWAS was analyzed by linear model, give model argument as 'linear'
GPC(OR = c(1.2,1.3,1.4,1.5,1.6), maf = c(0.1, 0.2, 0.3), N = 500, model ='linear')
```

```r
        GWAS Power Calculation       

_______________________________________ 
              Odds Ratio             
--------------------------------------- 
MAF|  1.2   1.3    1.4     1.5     1.6  
--------------------------------------- 
0.1| 0.01% 0.16%   1.3%   6.1%   18.38% 
0.2| 0.09% 1.81%  13.22% 42.78%  76.29% 
0.3| 0.26% 5.56%  32.45% 73.84%  95.54% 
--------------------------------------- 
```

#### 2. Logistics Model
```r
# If the GWAS was analyzed by logistics model, give model argument as 'binary'
GPC(OR = c(1.2,1.3,1.4,1.5,1.6), maf = c(0.1, 0.2, 0.3), N = 1500, model ='binary', Ncase = 500)
```

```r
       GWAS Power Calculation       

______________________________________ 
             Odds Ratio             
-------------------------------------- 
MAF|  1.2   1.3   1.4    1.5     1.6  
-------------------------------------- 
0.1|  0%   0.03%  0.22%  1.04%   3.51% 
0.2| 0.02% 0.31%  2.41%  10.32% 27.52% 
0.3| 0.05% 0.95%  7.08%  25.66% 54.37% 
-------------------------------------- 
```

