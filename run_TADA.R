# Copyright 2022 
# Original code and conception: Xin He
# Modifications and extensions: Jiebiao Wang, Lambertus Klei, and Jack Fu
# Distributed under terms of the MIT license.

library(GenomicRanges); library(stringr); library(openxlsx)
rm(list=ls())

### Change to working directory where the folder was cloned
source('functions.R')
### File name for supplementary table
supp_tab <- 'supplementary_table.xlsx'

########################################################################
### List of ACMG genes that were not released in certain cohorts
########################################################################
acmg <- c('ACTA2', 'ACTC1', 'APC', 'APOB', 'ATP7B', 'BMPR1A', 'BRCA1', 'BRCA2', 'CACNA1S', 'COL3A1', 'DSC2', 'DSG2', 'DSP', 'FBN1', 'GLA', 'KCNH2', 'KCNQ1', 'LDLR', 'LMNA', 'MEN1', 'MLH1', 'MSH2', 'MSH6', 'MUTYH', 'MYBPC3', 'MYH11', 'MYH7', 'MYL2', 'MYL3', 'NF2', 'OTC', 'PCSK9', 'PKP2', 'PMS2', 'PRKAG2', 'RB1', 'RET', 'RYR1', 'RYR2', 'SCN5A', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SMAD3', 'SMAD4', 'STK11', 'TGFBR1', 'TGFBR2', 'TMEM43', 'TNNI3', 'TNNT2', 'TP53', 'TPM1', 'VHL', 'WT1')

### Information on genes: constraint, priors
info <- read.xlsx(, sheet=8)

########################################################################
### Loading gene count tables
########################################################################
### SNV, denovo + inherited, ASC+SSC
dg_asc <- read.xlsx(supp_tab, sheet=5)

### SNV, denovo + inherited, SPARK
dg_spk <- read.xlsx(supp_tab, sheet=6)

### SNV, denovo, DDD
dg_ddd <- # See Data availability to access data

### SNV, case/control
dg_cc_dbs <- # See Data availability to access data
dg_cc_swe <- # See Data availability to access data
dg_cc <- read.xlsx(supp_tab, sheet=7)

########################################################################
### Loading CNV count tables
########################################################################
### CNV, denovo, ASC+SSC+SPARK
cnv_dn <- read.xlsx(supp_tab, sheet=9)
    cnv_dn <- cnv_dn[-which(cnv_dn$chr=="chr21" & (cnv_dn$end-cnv_dn$start)>30000000),] ### remove trisomy 21

### CNV, inherited, ASC+SSC+SPARK
cnv_in <- # See Data availability to access data

### CNV, case/control
cnv_cc <- # See Data availability to access data

########################################################################
### Setting sample sizes
########################################################################
## ASC B14 + SSC
n_prob_asc <- 7287; n_sib_asc <- 2348
## ASC B15-16
n_prob_asn <- 283; n_sib_asn <- 11
## Satterstrom et al lifted over
n_prob_lo <- 458; n_sib_lo <- 101 
## SPARK
n_prob_spk <- 6543; n_sib_spk <- 3032
## SPARK Pilot
n_prob_spp <- 465; n_sib_spp <- 0
## Case control
n1_case = 4863; n1_control = 5002 # DBS
n2_case <- 728; n2_control <- 3595 # SWE
## CNV
n_prob_cnv <- 13786; n_sib_cnv <- 5098
n_in_cnv <- 13182; n_in_cnv_sib <- 4866
n_case_cnv <- 684; n_control_cnv <- 10295

########################################################################
########################################################################
### Bayes Factor calculations
########################################################################
########################################################################
beta.dn <- 0.2

########################################################################
### SNV/indel
########################################################################

    ########################################################################
    ### DN PTV
    BF_dn_ptv_asc <- BF_DN_SNV(count_case=dg_asc$dn.ptv, count_con=dg_asc$dn.ptv.sib, n_case=n_prob_asc + n_prob_asn + n_prob_lo, n_con=n_sib_asc + n_sib_asn + n_sib_lo, mut=dg_asc$mut.ptv, gamma.dn=info$prior.dn.ptv, beta.dn=beta.dn)
    BF_dn_ptv_spk <- BF_DN_SNV(count_case=dg_spk$dn.ptv, count_con=dg_spk$dn.ptv.sib, n_case=n_prob_spk + n_prob_spp, n_con=n_sib_spk + n_sib_spp, mut=dg_spk$mut.ptv, gamma.dn=info$prior.dn.ptv, beta.dn=beta.dn)
    BF_dn_ptv <- BF_dn_ptv_asc*BF_dn_ptv_spk

    ### DN misB
    BF_dn_misB_asc <- BF_DN_SNV(count_case=dg_asc$dn.misb, count_con=dg_asc$dn.misb.sib, n_case=n_prob_asc + n_prob_asn + n_prob_lo, n_con=n_sib_asc + n_sib_asn + n_sib_lo, mut=dg_asc$mut.misb, gamma.dn=info$prior.dn.misb, beta.dn=beta.dn)
    BF_dn_misB_spk <- BF_DN_SNV(count_case=dg_spk$dn.misb, count_con=dg_spk$dn.misb.sib, n_case=n_prob_spk + n_prob_spp, n_con=n_sib_spk + n_sib_spp, mut=dg_spk$mut.misb, gamma.dn=info$prior.dn.misb, beta.dn=beta.dn)
    BF_dn_misB <- BF_dn_misB_asc*BF_dn_misB_spk

    ### DN misA
    BF_dn_misA_asc <- BF_DN_SNV(count_case=dg_asc$dn.misa, count_con=dg_asc$dn.misa.sib, n_case=n_prob_asc + n_prob_asn + n_prob_lo, n_con=n_sib_asc + n_sib_asc + n_sib_lo, mut=dg_asc$mut.misa, gamma.dn=info$prior.dn.misa, beta.dn=beta.dn)
    BF_dn_misA_spk <- BF_DN_SNV(count_case=dg_spk$dn.misa, count_con=dg_spk$dn.misa.sib, n_case=n_prob_spk + n_prob_spp, n_con=n_sib_spk + n_sib_spp, mut=dg_spk$mut.misa, gamma.dn=info$prior.dn.misa, beta.dn=beta.dn)
    BF_dn_misA <- BF_dn_misA_asc*BF_dn_misA_spk

    ########################################################################
    ### IN PTV
    BF_in_ptv <- BF_CC_SNV(count_case=dg_asc$in.t.ptv+dg_spk$in.t.ptv, count_con=dg_asc$in.u.ptv+dg_spk$in.u.ptv, n_case= n_prob_asc + n_prob_asn + n_prob_spk + n_prob_spp, n_con=n_prob_asc + n_prob_asn + n_prob_spk + n_prob_spp, mut=dg_asc$mut.ptv, gamma.cc=info$prior.in.ptv)
    ### IN misB
    BF_in_misB <- BF_CC_SNV(count_case=dg_asc$in.t.misb+dg_spk$in.t.misb, count_con=dg_asc$in.u.misb+dg_spk$in.u.misb,  n_case= n_prob_asc + n_prob_asn + n_prob_spk + n_prob_spp, n_con=n_prob_asc + n_prob_asn + n_prob_spk + n_prob_spp, mut=dg_asc$mut.misb, gamma.cc=info$prior.in.misb)
    ### IN misA
    BF_in_misA <- BF_CC_SNV(count_case=dg_asc$in.t.misa+dg_spk$in.t.misa, count_con=dg_asc$in.u.misa+dg_spk$in.u.misa,  n_case= n_prob_asc + n_prob_asn + n_prob_spk + n_prob_spp, n_con=n_prob_asc + n_prob_asn + n_prob_spk + n_prob_spp, mut=dg_asc$mut.misa, gamma.cc=info$prior.in.misa)

    ########################################################################
    ### CC PTV
    BF_cc_ptv <- BF_CC_SNV(count_case=dg_cc$case.ptv, count_con=dg_cc$control.ptv, n_case=n1_case+n2_case, n_con=n1_control+n2_control, mut=dg_asc$mut.ptv, gamma.cc=info$prior.cc.ptv)
    ### CC misB
    BF_cc_misB <- BF_CC_SNV(count_case=dg_cc$case.misb, count_con=dg_cc$control.misb, n_case=n1_case+n2_case, n_con=n1_control+n2_control, mut=dg_asc$mut.misb, gamma.cc=info$prior.cc.misb)
    ### CC misA
    BF_cc_misA <- BF_CC_SNV(count_case=dg_cc$case.misa, count_con=dg_cc$control.misa, n_case=n1_case+n2_case, n_con=n1_control+n2_control, mut=dg_asc$mut.misa, gamma.cc=info$prior.cc.misa)

########################################################################
### CNV
########################################################################
# size range, in number of constrained genes, to consider
cnv_size_range <- 1:8 

    ########################################################################
    ## DN
    cnv <- processCNV(cnv_dn, loeuf_threshold=0.6,  info=info)
    cnv_use <- cnv[which((is.na(cnv$gd_loci) | cnv$nahr==FALSE) & cnv$non_diploid_freq<.01 & cnv$num_genes %in% cnv_size_range & !is.na(cnv$Affected_Status)),]
        del_dup_adj <- table(cnv[,c("Affected_Status", "call")])
        del_dup_adj <- (del_dup_adj[2,1]/del_dup_adj[1,1])/(del_dup_adj[2,2]/del_dup_adj[1,2])

    del_use <- cnv_use[cnv_use$call=="DEL",]
    mut.pred.del <- sapply(cnv_size_range, function(x){ max(1, sum(del_use$num_genes==x &del_use$Affected_Status==1))/n_sib_cnv/(length(info$prior.dn.ptv)-x)})
    del_use$mut <- mut.pred.del[del_use$num_genes]

    BF_dn_del_prob <- BF_DN_CNV(cnv_size_range=cnv_size_range, cnv_use=del_use[del_use$Affected_Status==2,], info=info, prior=info$prior.dn.ptv, beta.dn=beta.dn, n=n_prob_cnv)
    BF_dn_del_sib <- BF_DN_CNV(cnv_size_range=cnv_size_range, cnv_use=del_use[del_use$Affected_Status==1,], info=info, prior=info$prior.dn.ptv, beta.dn=beta.dn, n=n_sib_cnv)
    BF_dn_del <- pmax(1, BF_dn_del_prob/BF_dn_del_sib)

    dup_use <- cnv_use[cnv_use$call=="DUP",]
    mut.pred.dup <- sapply(cnv_size_range, function(x){ max(1, sum(dup_use$num_genes==x &dup_use$Affected_Status==1))/n_sib_cnv/(length(info$prior.dn.ptv)-x)})
    dup_use$mut <- mut.pred.dup[dup_use$num_genes]

    BF_dn_dup_prob <- BF_DN_CNV(cnv_size_range=cnv_size_range, cnv_use=dup_use[dup_use$Affected_Status==2,], info=info, prior=info$prior.dn.ptv, beta.dn=beta.dn, n=n_prob_cnv, del_dup_adj=del_dup_adj)
    BF_dn_dup_sib <- BF_DN_CNV(cnv_size_range=cnv_size_range, cnv_use=dup_use[dup_use$Affected_Status==1,], info=info, prior=info$prior.dn.ptv, beta.dn=beta.dn, n=n_sib_cnv, del_dup_adj=del_dup_adj)
    BF_dn_dup <- pmax(1, BF_dn_dup_prob/BF_dn_dup_sib)

    ########################################################################
    ## CC
    cnv_size_range <- 1:8

    cnv <- processCNV(cnv_cc[which(cnv_cc$sf<.01 & (is.na(cnv_cc$gd_loci) | cnv_cc$nahr==FALSE)),], loeuf_threshold=0.6, info=info)
    cnv_use <- cnv[cnv$sf<.01 & cnv$num_genes %in% cnv_size_range,]

    BF_cc_del <- BF_CC_CNV(cnv_size_range, cnv_use=cnv_use[cnv_use$call=="DEL",], n_case=n_case_cnv, n_con=n_control_cnv, info=info, prior=info$prior.cc.ptv, del_dup_adj=1, nu=5000)
    BF_cc_dup <- BF_CC_CNV(cnv_size_range, cnv_use=cnv_use[cnv_use$call=="DUP",], n_case=n_case_cnv, n_con=n_control_cnv, info=info, prior=info$prior.cc.ptv, del_dup_adj=del_dup_adj, nu=5000)

    ########################################################################
    ## IN
    cnv <- processCNV(cnv_in[which(cnv_in$sf<.01 & (is.na(cnv_in$gd_loci) | cnv_in$nahr==FALSE) & cnv_in$Affected_Status_Child==2),], loeuf_threshold=0.6, info=info)
    cnv_use <- cnv[cnv$sf<.01 & cnv$num_genes %in% cnv_size_range & is.na(cnv$gd_loci),]

    BF_in_del <- BF_CC_CNV(cnv_size_range, cnv_use=cnv_use[cnv_use$call=="DEL",], n_case=n_in_cnv, n_con=n_in_cnv, info=info, prior=info$prior.in.ptv, del_dup_adj=1, nu=5000)
    BF_in_dup <- BF_CC_CNV(cnv_size_range, cnv_use=cnv_use[cnv_use$call=="DUP",], n_case=n_in_cnv, n_con=n_in_cnv, info=info, prior=info$prior.in.ptv, del_dup_adj=del_dup_adj, nu=5000)

########################################################################
### Accounting for genes with ACMG masking
########################################################################
acmg_ind <- which(info$gene %in% acmg)

    BF_dn_ptv_spk_acmg <- BF_DN_SNV(count_case=dg_spk$dn.ptv, count_con=dg_spk$dn.ptv.sib, n_case= n_prob_spp, n_con=0, mut=dg_spk$mut.ptv, gamma.dn=info$prior.dn.ptv, beta.dn=beta.dn)
        BF_dn_ptv_spk[acmg_ind] <- BF_dn_ptv_spk_acmg[acmg_ind]
    BF_dn_ptv <- BF_dn_ptv_asc*BF_dn_ptv_spk

    BF_dn_misB_spk_acmg <- BF_DN_SNV(count_case=dg_spk$dn.misb, count_con=dg_spk$dn.misb.sib,n_case= n_prob_spp, n_con=0, mut=dg_spk$mut.misb, gamma.dn=info$prior.dn.misb, beta.dn=beta.dn)
        BF_dn_misB_spk[acmg_ind] <- BF_dn_misB_spk_acmg[acmg_ind]
    BF_dn_misB <- BF_dn_misB_asc*BF_dn_misB_spk

    BF_dn_misA_spk_acmg <- BF_DN_SNV(count_case=dg_spk$dn.misa, count_con=dg_spk$dn.misa.sib,n_case= n_prob_spp, n_con=0, mut=dg_spk$mut.misa, gamma.dn=info$prior.dn.misa, beta.dn=beta.dn)
        BF_dn_misA_spk[acmg_ind] <- BF_dn_misA_spk_acmg[acmg_ind]
    BF_dn_misA <- BF_dn_misA_asc*BF_dn_misA_spk

    ### CC PTV
    BF_cc_ptv_acmg <- BF_CC_SNV(count_case=dg_cc_swe$case.ptv, count_con=dg_cc_swe$control.ptv, n_case=n2_case, n_con=n2_control, mut=dg_asc$mut.ptv, gamma.cc=info$prior.cc.ptv)
        BF_cc_ptv[acmg_ind] <- BF_cc_ptv_acmg[acmg_ind]

    ### CC misB
    BF_cc_misB_acmg <- BF_CC_SNV(count_case=dg_cc_swe$case.misb, count_con=dg_cc_swe$control.misb, n_case=n2_case, n_con=n2_control, mut=dg_asc$mut.misb, gamma.cc=info$prior.cc.misb)
        BF_cc_misB[acmg_ind] <- BF_cc_misB_acmg[acmg_ind]

    ### CC misA
    BF_cc_misA_acmg <- BF_CC_SNV(count_case=dg_cc_swe$case.misa, count_con=dg_cc_swe$control.misa, n_case=n2_case, n_con=n2_control, mut=dg_asc$mut.misa, gamma.cc=info$prior.cc.misa)
        BF_cc_misA[acmg_ind] <- BF_cc_misA_acmg[acmg_ind]
    
########################################################################
### Integrating DDD - SNV/indel
########################################################################

### DN PTV
BF_dn_ptv_ddd <- BF_DN_SNV(count_case=dg_ddd$dn.ptv, count_con=0, n_case=31058, n_con=0, mut=dg_asc$mut.ptv, gamma.dn=info$prior.dn.ptv, beta.dn=beta.dn)

### DN misB
BF_dn_misB_ddd <- BF_DN_SNV(count_case=dg_ddd$dn.misb, count_con=0, n_case=31058, n_con=0, mut=dg_asc$mut.misb, gamma.dn=info$prior.dn.misb, beta.dn=beta.dn)

### DN misA
BF_dn_misA_ddd <- BF_DN_SNV(count_case=dg_ddd$dn.misa, count_con=0, n_case=31058, n_con=0, mut=dg_asc$mut.misa, gamma.dn=info$prior.dn.misa, beta.dn=beta.dn)


########################################################################
########################################################################
### TADA outcome calculation
########################################################################
########################################################################

### Aggregate ASD BF
BF_asd <- cbind(pmax(BF_dn_ptv_asc*BF_dn_ptv_spk*BF_in_ptv*BF_cc_ptv, 1),
    pmax(BF_dn_misB_asc*BF_dn_misB_spk*BF_in_misB*BF_cc_misB, 1),
    pmax(BF_dn_misA_asc*BF_dn_misA_spk*BF_in_misA*BF_cc_misA, 1),
    pmax(BF_dn_del*BF_cc_del*BF_in_del, 1), pmax(BF_dn_dup*BF_cc_dup*BF_in_dup, 1))
# Only use del/dup if there is snv/indel evidence (Tot BF > 2) 
BF_asd[apply(BF_asd[,1:3], 1, max)<5,4:5] <- 1

### Aggregate DD BF
BF_ddd <- cbind(pmax(1, BF_dn_ptv_ddd), pmax(BF_dn_misB_ddd, 1), pmax(BF_dn_misA_ddd, 1))

### Aggregate NDD BF
BF_asd_ddd <- cbind(BF_asd, BF_ddd)

### FDR calculation
qval_asd <- Bayesian.FDR(apply(BF_asd, 1, prod), pi0 = 1 - 0.05)
    names(qval_asd) <- info$gene
    ### TADA-ASD significant genes at various thresholds
    sapply(c(0.001, 0.01, 0.05, 0.1), function(x) sum(qval_asd<x))

qval_ddd <- Bayesian.FDR(apply(BF_ddd, 1, prod), pi0 = 1 - 0.05)
    names(qval_ddd) <- info$gene
    ### TADA-DD significant genes at various thresholds
    sapply(c(0.001, 0.01, 0.05, 0.1), function(x) sum(qval_ddd<x))

qval_asd_ddd <- Bayesian.FDR(apply(BF_asd_ddd, 1, prod), pi0=1-.05)
    names(qval_asd_ddd) <- info$gene
    ### TADA-NDD significant genes at various thresholds
    sapply(c(0.001, 0.01, 0.05, 0.1), function(x) sum(qval_asd_ddd<x))

### back-transform FDR to p-values
ec_threshold <- .05/length(qval_asd)
pval_asd <- q2p(qval_asd)
    sum(pval_asd <= ec_threshold)
pval_ddd <- q2p(qval_ddd)
    sum(pval_ddd <= ec_threshold)
pval_asd_ddd <- q2p(qval_asd_ddd)
    sum(pval_asd_ddd <= ec_threshold)
