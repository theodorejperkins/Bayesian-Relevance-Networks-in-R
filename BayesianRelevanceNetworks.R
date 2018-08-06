##############################################################################
# This R source code is in support of the publication "Uncovering Robust 
# Patterns of MicroRNA Co-Expression across Cancers using Bayesian Relevance 
# Networks" by Ramachandran, Sanchez-Taltavull, and Perkins, PLoS ONE, 
# Vol. 12, No. 8, Art. e0183103 (also appeared at GLBIO 2017, where it won 
# "Outstanding Presentation" prize).
# 
# Questions or comments: theodore.j.perkins@gmail.com  or rpara26@gmail.com  
##############################################################################

#######################################################
# FUNCTION BayesianCorrelation_Grouped
# #####################################################
# 
# This function computes the grouped Bayesian correlations between all pairs of
# m entities across n conditions. The m-by-n "ReadCounts" input matrix specifies
# the numbers of reads for each entity (rows) and condition (columns). The 
# second input, "Groups", is a 1-by-n vector of group numbers, specifying to 
# which group each condition belongs. For instance, if the first two conditions 
# are group 1, second three conditions are group 3, and third three conditions 
# are group 2, we would have Groups = [1 1 3 3 3 2 2 2]. The optional 3rd 
# argument, if true, tells the function not to compute the second half of the 
# covariance term (covariance of uncertainties in levels, averaged across 
# conditions & groups). The term tends to be very small compared to everything 
# else, yet slow to compute. Further, it is unchanged in permutation 
# computations, so in most cases, it can be safely ignored. The final result is 
# returned in the m-by-m matrix Corrs.

BayesianCorrelation_Grouped <- function(ReadCounts,Groups,SkipCovU){

if (missing(SkipCovU)) SkipCovU = 0

# How big are things?
m <- nrow(ReadCounts)
s <- ncol(ReadCounts)

UniqGroups = unique(Groups)
g = length(UniqGroups)

# Priors
PriorAlphas_ms <- matrix(rep(rep(1,s),m),nrow=m,byrow=TRUE)/(m-1)

# Compute posteriors and concentration parameters (sums of alphas across entities)
PosteriorAlphas_ms = ReadCounts + PriorAlphas_ms
TotalAlphas_s = colSums(PosteriorAlphas_ms)
TotalAlphas_ms = matrix(rep(TotalAlphas_s,m), nrow=m, byrow=TRUE)

# Posterior means and variances by sample
MeanU_ms = PosteriorAlphas_ms/TotalAlphas_ms
VarU_ms = PosteriorAlphas_ms*(TotalAlphas_ms-PosteriorAlphas_ms)/(TotalAlphas_ms*TotalAlphas_ms*(TotalAlphas_ms+1))

# Posterior means and variances by group
MeanU_mg = NULL
VarU_mg = NULL
for (i in 1:g){
I = which(Groups==UniqGroups[i])
tmprowMeans <- rowMeans(MeanU_ms[,I])
MeanU_mg = cbind(MeanU_mg, tmprowMeans)
X = MeanU_ms[,I]-matrix(rep(tmprowMeans,length(I)), ncol=length(I))
VarU_mg = cbind(VarU_mg, (rowMeans(VarU_ms[,I])+rowMeans(X^2))/length(I))
}

# Variance across groups
MeanGMeanU_m = rowMeans(MeanU_mg)
X = MeanU_mg-matrix(rep(MeanGMeanU_m,g), ncol=g)
VarGMeanU_m = rowMeans(X^2)

# Total variance
MeanGVarU_m = rowMeans(VarU_mg)
VarGU_m = VarGMeanU_m + MeanGVarU_m

# Mean across groups of covariance within each group
MeanGCovU_mm = mat.or.vec(m,m)
if (!SkipCovU){
        for (i in 1:g){
                TempCov = mat.or.vec(m,m)
                J = which(Groups==UniqGroups[i])
                for (j in 1:length(J)){
                        s = J[j]
                        TempCov = TempCov - MeanU_ms[,s]%*%t(MeanU_ms[,s])/(TotalAlphas_s[s]+1)
                }
                TempCov = TempCov/length(J)
                MeanGCovU_mm = MeanGCovU_mm + TempCov
        }
        MeanGCovU_mm = MeanGCovU_mm/g
}

# Covariance across groups of mean within each group
CovGMeanU_mm = X%*%t(X)/g

# Total covariance
CovGU_mm = CovGMeanU_mm + MeanGCovU_mm

# Correlation
S_m = 1/sqrt(VarGU_m)
SS_mm = S_m%*%t(S_m)
Corrs = CovGU_mm * SS_mm
Corrs[seq(1, m^2, by=m+1)] = 1

return(Corrs)
}
#--------------------------------------------------
# End of FUNCTION BayesianCorrelation_Grouped
#--------------------------------------------------

#######################################################
# FUNCTION BayesianPermutation_Grouped
#######################################################
# 
# This function estimates a null distribution for corrleations computed by the
# function BayesianCorrelation_Grouped. The inputs are ReadCounts (m-by-n
# matrix), Groups (1-by-n vector), Repeats, a positive integer, and SkipCovU.
# The ReadCounts, Groups and SkipCovU inputs have the same meaning as for the 
# function BayesianCorrelation_Grouped. The Repeats input specifies the number
# of random permutations to test. The output is an m-by-m-by-Repeats matrix of
# correlations computed from permutations. These can be used as estimates of
# null distributions for each pair of entities, or can be combined to form a
# single, overall null distribution.

BayesianPermutation_Grouped <- function(ReadCounts, Groups, Repeats, SkipCovU){
        
        print(" ",quote=FALSE)
        print("------------------",quote=FALSE)
        print("Bayesian Permutation ... ",quote=FALSE)
        print("------------------",quote=FALSE)
        print(" ",quote=FALSE)
        
        if (missing(SkipCovU)) SkipCovU = 0
        
        # How big are things?
        m <- nrow(ReadCounts)
        s <- ncol(ReadCounts)
        
        UniqGroups = unique(Groups)
        g = length(UniqGroups)
        
        # Priors
        PriorAlphas_ms <- matrix(rep(rep(1,s),m),nrow=m,byrow=TRUE)/(m-1)
        
        # Compute posteriors and concentration parameters (sums of alphas across entities)
        PosteriorAlphas_ms = ReadCounts + PriorAlphas_ms
        TotalAlphas_s = colSums(PosteriorAlphas_ms)
        TotalAlphas_ms = matrix(rep(TotalAlphas_s,m), nrow=m, byrow=TRUE)
        
        # Posterior means and variances by sample
        MeanU_ms = PosteriorAlphas_ms/TotalAlphas_ms
        VarU_ms = PosteriorAlphas_ms*(TotalAlphas_ms-PosteriorAlphas_ms)/(TotalAlphas_ms*TotalAlphas_ms*(TotalAlphas_ms+1))
        
        # Posterior means and variances by group
        MeanU_mg = NULL
        VarU_mg = NULL
        for (i in 1:g){
                I = which(Groups==UniqGroups[i])
                tmprowMeans <- rowMeans(MeanU_ms[,I])
                MeanU_mg = cbind(MeanU_mg, tmprowMeans)
                X = MeanU_ms[,I]-matrix(rep(tmprowMeans,length(I)), ncol=length(I))
                VarU_mg = cbind(VarU_mg, (rowMeans(VarU_ms[,I])+rowMeans(X^2))/length(I))
        }
        
        # Variance across groups
        MeanGMeanU_m = rowMeans(MeanU_mg)
        X = MeanU_mg-matrix(rep(MeanGMeanU_m,g), ncol=g)
        VarGMeanU_m = rowMeans(X^2)
        
        # Total variance
        MeanGVarU_m = rowMeans(VarU_mg)
        VarGU_m = VarGMeanU_m + MeanGVarU_m
        
        # Mean across groups of covariance within each group
        MeanGCovU_mm = mat.or.vec(m,m)
        if (!SkipCovU){
                for (i in 1:g){
                        TempCov = mat.or.vec(m,m)
                        J = which(Groups==UniqGroups[i])
                        for (j in 1:length(J)){
                                s = J[j]
                                TempCov = TempCov - MeanU_ms[,s]%*%t(MeanU_ms[,s])/(TotalAlphas_s[s]+1)
                        }
                        TempCov = TempCov/length(J)
                        MeanGCovU_mm = MeanGCovU_mm + TempCov
                }
                MeanGCovU_mm = MeanGCovU_mm/g
        }
        
        # Permuted covariance and correlation
        S_m = 1/sqrt(VarGU_m)
        SS_mm = S_m%*%t(S_m)
        Corrs_mmr = array(0, dim=c(m,m,Repeats))
        for (r in 1:Repeats){
                print(paste("Repeat ", r, sep=""))
                # Permute one copy of X
                TempX = X
                for (i in 1:m){
                        TempX[i,] = TempX[i,sample(g)]
                }
                
                # Permuted covariance
                CovGMeanU_mm = X%*%t(TempX)/g
                CovGU_mm = CovGMeanU_mm + MeanGCovU_mm
                # Permuted correlation
                Corrs = CovGU_mm * SS_mm
                Corrs[seq(1, m^2, by=m+1)] = 1
                Corrs_mmr[,,r] = Corrs
        }
        return(Corrs_mmr)
}
#--------------------------------------------------
# End of FUNCTION BayesianPermutation_Grouped
#--------------------------------------------------

#######################################################
# FUNCTION PearsonCorrelation_Grouped
#######################################################

PearsonCorrelation_Grouped <- function(ReadCounts, Groups){
        
# This function computes the grouped Pearson correlations between all pairs of m
# entities across n conditions. The m-by-n ReadCounts input matrix specifies the
# numbers of reads for each entity (rows) and condition (columns). The second
# input, Groups, is a 1-by-n vector of group numbers, specifying to which group 
# each condition belongs. For instance, if the first two conditions are group 1,
# second three conditions are group 3, and third three conditions are group 2,
# we would have Groups = [1 1 3 3 3 2 2 2]. The correlation is computed by
# normalizing read counts by the total reads in each column, then averaging
# those within groups, and computing the correlation across groups. The answer
# is return in the m-by-m matrix Corrs.
        
        # Sizes of things
        m=nrow(ReadCounts)
        n=ncol(ReadCounts)
        
        UniqGroups = unique(Groups)
        g = length(UniqGroups)
        
        # Empirical fractions
        F_mn = mat.or.vec(m,n)
        for (j in 1:n){
                F_mn[,j] = ReadCounts[,j]/sum(ReadCounts[,j])
        }
        
        # Group means
        MeanS_mg = mat.or.vec(m,g)
        for (i in 1:g){
                I = which(Groups==UniqGroups[i])
                MeanS_mg[,i] = rowMeans(F_mn[,I])
        }
        
        # Entity means across groups
        MeanGS_m = rowMeans(MeanS_mg)
        
        # Entity variances across groups
        X_mg = MeanS_mg - matrix(rep(MeanGS_m,g), ncol=g)    
        VarG_m = rowMeans(X_mg^2)
        OneByStd_m = 1/sqrt(VarG_m)
        OneByStd_mm = OneByStd_m %*% t(OneByStd_m)
        
        # Entity covariances across groups
        CovG_mm = X_mg %*% t(X_mg)/g
        
        # Correlations
        Corrs = CovG_mm*OneByStd_mm
        Corrs[seq(1, m^2, by=m+1)] = 1
        Corrs[which(is.nan(Corrs))]=0
        return(Corrs)
}
#--------------------------------------------------
# End of FUNCTION PearsonCorrelation_Grouped
#--------------------------------------------------

#######################################################
# FUNCTION PearsonPermutation_Grouped
#######################################################

PearsonPermutation_Grouped <- function(ReadCounts, Groups, Repeats){
        
# This function estimates a null distribution for corrleations computed by the 
# function PearsonCorrelation_Grouped. The inputs are ReadCounts (m-by-n 
# matrix), Groups (1-by-n vector), and repeats, a positive integer. The 
# ReadCounts and Groups inputs have the same meaning as for the function 
# PearsonCorrelation_Grouped. The final input specifies the number of random 
# permutations to test. The output is an m-by-m-by-Repeats matrix of 
# correlations computed from permutations. These can be used as estimates of 
# null distributions for each pair of entities, or can be combined to form a 
# single, overall null distribution.
        
        print(" ",quote=FALSE)
        print("------------------",quote=FALSE)
        print("Pearson Permutation ... ",quote=FALSE)
        print("------------------",quote=FALSE)
        print(" ",quote=FALSE)
        
        # Sizes of things
        m=nrow(ReadCounts)
        n=ncol(ReadCounts)
        
        UniqGroups = unique(Groups)
        g = length(UniqGroups)
        
        # Empirical fractions
        F_mn = mat.or.vec(m,n)
        for (j in 1:n){
                F_mn[,j] = ReadCounts[,j]/sum(ReadCounts[,j])
        }
        
        # Group means
        MeanS_mg = mat.or.vec(m,g)
        for (i in 1:g){
                I = which(Groups==UniqGroups[i])
                MeanS_mg[,i] = rowMeans(F_mn[,I])
        }
        
        # Entity means across groups
        MeanGS_m = rowMeans(MeanS_mg)
        
        # Entity variances across groups
        X_mg = MeanS_mg - matrix(rep(MeanGS_m,g), ncol=g)    
        VarG_m = rowMeans(X_mg^2)
        OneByStd_m = 1/sqrt(VarG_m)
        OneByStd_mm = OneByStd_m %*% t(OneByStd_m)
        
        # Permuted covariances and correlations
        Corrs_mmr = array(0, dim=c(m,m,Repeats))
        for (r in 1:Repeats){
                print(paste("Repeat ", r, sep=""))
                # Permute X temporarily
                TempX_mg = X_mg
                for (i in 1:m){
                        TempX_mg[i,] = TempX_mg[i,sample(g)]
                }
                # Permuted entity covariances across groups
                TempCov = X_mg %*% t(TempX_mg)/g
                # Correlations
                TempCorrs = TempCov * OneByStd_mm
                TempCorrs[seq(1, m^2, by=m+1)] = 1
                Corrs_mmr[,,r] = TempCorrs
        }
        Corrs_mmr[which(is.nan(Corrs_mmr))]=0
        return(Corrs_mmr)
}
#--------------------------------------------------
# End of FUNCTION PearsonPermutation_Grouped
#--------------------------------------------------

#######################################################
# FUNCTION FDRAnalysis
#######################################################

FDRAnalysis <- function(Corrs, PermCorrs){
        
        # This function performs a false discovery rate analysis of m-by-m correlation
        # matrix Corrs, in comparison with the permutation-based m-by-m-by-Repeats
        # correlation matrix PermCorrs. It computes how many above-diagonal entries of
        # Corrs are above different possible correlation thresholds (namely -1:0.01:1).
        # It looks at the empirical fraction of permuted correlations above each of
        # those thresholds, on the basis of which it computes expected numbers of false
        # positives. And from that, it estimates false discovery rate. These information
        # provide guidance to the user for selecting a correlation cutoff to form a
        # relevance network.
        
        m <- dim(Corrs)[1]
        # dummy1 <- dim(Corrs)[2]
        #
        # dummy2 <- dim(PermCorrs)[1]
        # dummy3 <- dim(PermCorrs)[2]
        Repeats <- dim(PermCorrs)[3]
        
        # What is the set of correlation thresholds that we will test?
        FDR.CorrThresh = seq(-1, 1, by=0.01)
        FDR.NCorrAbove <- NULL
        FDR.FPermAbove <- NULL
        FDR.EFalsePos <- NULL
        FDR.EFDR <- NULL
        
        # For each threshold...
        for (i in 1:length(FDR.CorrThresh)){
                # Threshold
                CT = FDR.CorrThresh[i]
                # How many above-diagonal entries are above that threshold?
                FDR.NCorrAbove[i] = (sum(Corrs>=CT)-m)/2
                # What fraction of permutations are above that threshold?
                FDR.FPermAbove[i] = (sum(PermCorrs>=CT)-m*Repeats)/(Repeats*m*(m-1))
                # Expected false positives
                FDR.EFalsePos[i] = FDR.FPermAbove[i]*m*(m-1)/2
                # Estimated false discovery rate
                FDR.EFDR[i] = FDR.EFalsePos[i]/FDR.NCorrAbove[i]
        }
        FDR <- data.frame(CorrThresh=FDR.CorrThresh, NCorrAbove=FDR.NCorrAbove, FPermAbove=FDR.FPermAbove, EFalsePos=FDR.EFalsePos, EFDR=FDR.EFDR)
        FDR$EFDR[FDR$EFDR>1] <- 1
        return(FDR)
}
#--------------------------------------------------
# End of FUNCTION FDRAnalysis
#--------------------------------------------------

#############################################################
# T E S T *** S C R I P T ***
#############################################################
ReadCounts = matrix(c(10,13,12,30,35,25,19,20,22,100,110,106,
                      400,300,350,90,100,92,16,14,12,33,40,35,
                      100,110,120,300,333,290,200,225,212,800,810,820,
                      900,920,919,200,210,211,50,60,55,65,70,75,
                      0,0,0,0,0,0,0,0,0,0,0,1,
                      0,0,0,0,0,0,0,0,0,0,0,2,
                      1,2,3,4,5,6,7,8,9,10,11,12), nrow=7, byrow=TRUE)

Groups = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4)

# Bayesian correlations
BCorr = BayesianCorrelation_Grouped(ReadCounts,Groups)

# Permuted Bayesian correlations
BPerm = BayesianPermutation_Grouped(ReadCounts,Groups,100)

# False discovery rate analysis
BFDR = FDRAnalysis(BCorr,BPerm)

# Suppose we choose a correlation threshold of 0.8, the following would be
# the links in the Bayesian Relevance Network.
BI <- which(BCorr>=0.8, arr.ind=TRUE)[,1]
BJ <- which(BCorr>=0.8, arr.ind=TRUE)[,2]
K = which(BJ > BI)
BI = BI[K]
BJ = BJ[K]

# Pearson correlations
PCorr = PearsonCorrelation_Grouped(ReadCounts,Groups)

# Permuted Pearson correlations
PPerm = PearsonPermutation_Grouped(ReadCounts,Groups,100)

# False discovery rate analysis
PFDR = FDRAnalysis(PCorr,PPerm)

# Suppose we choose a correlation threshold of 0.8, the following would be
# the links in the Pearson Relevance Network.
PI <- which(PCorr>=0.8, arr.ind=TRUE)[,1]
PJ <- which(PCorr>=0.8, arr.ind=TRUE)[,2]
K = which(PJ > PI)
PI = PI[K]
PJ = PJ[K]

