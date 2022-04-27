# Copyright 2022 
# Original code and conception: Xin He
# Modifications and extensions: Jiebiao Wang, Lambertus Klei, and Jack Fu
# Distributed under terms of the MIT license.

## Functions used to apply the TADA framework

################################################################################################
### log Bayes factor calculation for de novo variant contribution to a gene
## Input:
# x.dn: the de novo variant count 
# n.dn: the sample size (number of families)
# mu: the mutation rate (of this type of mutational events in this gene)
# Prior distribution of risk: gamma ~ Gamma(gamma.dn*beta.dn, beta.dn)
## Output:
# lBF: log Bayes factor evidence for de novo variants of this type in this gene
################################################################################################
log.bayes.factor.dn <- function(x.dn, n.dn, mu,gamma.dn, beta.dn){
  marg.lik0 <- dpois(x.dn, 2*n.dn*mu, log=TRUE)
  marg.lik1 <- dnbinom(x.dn, gamma.dn*beta.dn, beta.dn/(beta.dn+2*n.dn*mu), log=TRUE)
  lBF <- marg.lik1-marg.lik0
  return (lBF=lBF)
}

################################################################################################
### Bayes factor calculation for case-control/inherited contribution to a gene
## Input: 
# x.cc is a vector of case/transmitted and control/untransmitted variant counts in this gene
# n.cc is the number of case and control samples (or number of probands when evaluating inherited variants)
# gamma.cc is a vector of the prior mean on risk
# Prior distribution of q|H1: Gamma(rho1, nu1)
# Prior distribution of q|H0: Gamma(rho0, nu0)
## Output:
# BF: Bayes factor evidence for case-control/inherited contribution to a gene
################################################################################################
bayes.factor.cc <- function(x.cc, n.cc, gamma.cc, rho1, nu1, rho0, nu0){
  marglik0.cc <- evidence.null.cc(x.cc, n.cc, rho0, nu0)
  marglik1.cc <- evidence.alt.cc(x.cc, n.cc, gamma.cc, rho1, nu1)
  BF.cn <- marglik1.cc$cn / marglik0.cc$cn
    ## Assuming null and alt of the control model is identical
    if(is.na(BF.cn)){BF.cn<-1} 
  BF.ca <- marglik1.cc$ca / marglik0.cc$ca
  BF <- BF.cn * BF.ca
  return(BF=BF)
}

################################################################################################
### Helper function to bayes.factor.cc
### Calculates the evidence of case-control/inherited counts under null model
## Input: 
# x.cc is a vector of case/transmitted and control/untransmitted variant counts in this gene
# n.cc is the number of case and control samples (or number of probands when evaluating inherited variants)
# Prior distribution of q|H0: Gamma(rho0, nu0)
## Output:
# cn: marginal likelihood of control/untransmitted variant counts under null model
# ca: marginal likelihood of case/transmitted variant counts under null model
################################################################################################
evidence.null.cc <- function(x.cc, n.cc, rho0, nu0) {
  marglik0.ctrl.log <- log(dnbinom(x.cc$cn, rho0, nu0/(nu0+n.cc$cn)))
  marglik0.case.log <- log(dnbinom(x.cc$ca, rho0+x.cc$cn, (nu0+n.cc$cn)/(nu0+n.cc$cn+n.cc$ca)))
  marglik0.log <- marglik0.ctrl.log + marglik0.case.log
  
  return (list(cn=exp(marglik0.ctrl.log), ca=exp(marglik0.case.log)))#, total=exp(marglik0.log)))
}

################################################################################################
### Helper function to bayes.factor.cc
### Calculates the evidence of case-control/inherited counts under alternate model
## Input: 
# x.cc is a vector of case/transmitted and control/untransmitted variant counts in this gene
# n.cc is the number of case and control samples (or number of probands when evaluating inherited variants)
# gamma.cc is a vector of the prior mean on risk
## Output:
# cn: marginal likelihood of control/untransmitted variant counts under alternate model
# ca: marginal likelihood of case/transmitted variant counts under alternate model
################################################################################################
evidence.alt.cc <- function(x.cc, n.cc, gamma.cc, rho1, nu1){
  marglik1.ctrl <- dnbinom(x.cc$cn, rho1, nu1/(nu1+n.cc$cn))
  marglik1.case <- dnbinom(x.cc$ca, rho1+x.cc$cn, (nu1+n.cc$cn)/(nu1+n.cc$ca*gamma.cc+n.cc$cn))
  marglik1 <- marglik1.ctrl * marglik1.case
  return (list(cn=marglik1.ctrl, ca=marglik1.case))#, total=marglik1))
}

################################################################################################
### Wrapper function to obtain de novo Bayes Factors for SNV/indels, over a list of genes for a given variant type
### Makes calls to log.bayes.factor.dn for each gene
## Input: 
# count_case: a vector of de novo variant counts over each gene in cases
# count_con: a vector of de novo variant counts over each gene in controls
# n_case: the number of cases
# n_con: the number of controls
# mut: a vector of mutation rates for each gene
# Prior distribution of RR: gamma ~ Gamma(gamma.dn*beta.dn, beta.dn)
## Output:
# BF: a vector of Bayes factors of de novo evidence from SNV/indel for each gene for a given variant type
################################################################################################
BF_DN_SNV <- function(count_case, count_con, n_case, n_con, mut, gamma.dn, beta.dn=.2){
  ### If there are siblings, evaluate BF evidence for siblings to offset evidence for probands
  if(n_con>0){
    BF_con <- exp(sapply(1:length(mut), function(i) log.bayes.factor.dn(x.dn=count_con[i], n.dn=n_con, mu=mut[i], beta.dn=beta.dn, gamma.dn=gamma.dn[i])))
  }else{
    BF_con <- rep(1, length(mut))
  }
  BF_con[BF_con<1] <- 1

  ### BF evidence calculation for probands
  BF_case <- exp(sapply(1:length(mut), function(i) log.bayes.factor.dn(x.dn=count_case[i], n.dn=n_case, mu=mut[i], beta.dn=beta.dn, gamma.dn=gamma.dn[i])))
    BF_case[BF_case<1] <- 1
  
  ### Combine proband and sibling BF evidence
  ### Enforce BF back to 1 if no observations in probands or siblings
  BF <- pmax(BF_case/BF_con, 1) 
    BF[which(count_case==0 & count_con==0)] <- 1
    BF[is.na(BF)] <- 1
  return(BF)
}

################################################################################################
### Wrapper function to obtain case-control/inherited Bayes Factors for SNV/indels, over a list of genes for a given variant type
### Makes calls to bayes.factor.cc for each gene
## Input: 
# count_case: a vector of case/transmitted variant counts over each gene in cases
# count_con: a vector of control/untransmitted variant counts over each gene in controls
# n_case: the number of cases
# n_con: the number of controls (case-control analysis) or number of cases (inherited analysis)
# mut: a vector of mutation rates for each gene
# gamma.cc: a vector of the prior mean on risk
# nu: a hyperparameter that controls nu0 and nu1 in determining prior distributions, passed to bayes.factor.cc
## Output:
# BF: a vector of Bayes factors of case-control/inherited evidence from SNV/indel for each gene for a given variant type
################################################################################################
BF_CC_SNV <- function(count_case, count_con, n_case, n_con, mut, gamma.cc, nu=5000){
  rho.in = nu * sum(count_con, na.rm=T) / (2*length(mut)*n_con)
  rho.in = as.numeric(rho.in)*mut/mean(mut, na.rm=TRUE)
  BF = sapply((1:length(mut)), function(i){
    bayes.factor.cc(x.cc = data.frame(ca=count_case[i], cn=count_con[i]), n.cc=data.frame(ca=n_case, cn=n_con), gamma.cc=gamma.cc[i], rho1=rho.in[i], nu1=nu, rho0=rho.in[i], nu0=nu)
    })

  ### Enforce BF back to 1 if no mutation rate or no observed variant counts
  BF[mut==0] <- 1
  BF[which((count_case+count_con)==0)] <- 1
  return(BF) 
}

################################################################################################
### Bayes factor calculation for de novo CNV contribution to a gene
## Input: 
# cnv_size_range: a vector of size range of CNV events to considered (# of constrained genes)
# cnv_use:  the set of cnvs to be fed into the model (for format refer to Supplementary Table 9 and run_TADA.R)
# info: matrix that contains gene level information (Supplementary Table 8: gene names, constraint)
# prior: a vector of the prior risk for each gene
# Prior distribution of risk: gamma ~ Gamma(prior*beta.dn, beta.dn)
# n: the number of samples
# del_dup_adj: the adjustment for duplications being in general less penetrant than deletions
## Output:
# BF: a vector of Bayes factors of de novo evidence from CNVs for each gene for a given variant type
################################################################################################
BF_DN_CNV <- function(cnv_size_range, cnv_use, info, prior, beta.dn=0.2, n, del_dup_adj=1){
    raw_lBF <- do.call(rbind, lapply(cnv_size_range, function(i, cnv_use, n, del_dup_adj, prior){
        ### Subset to cnv of a gene-size
        cnv_sub <- cnv_use[cnv_use$num_genes==i, ]
        if(dim(cnv_sub)[1]>0){
            cnv_sub$relative_risk <- sapply(cnv_sub$genes, function(x) (sum(prior[match(x, info$gene_gencodeV33)])-1)/del_dup_adj+1)

            genes_combo <- table(cnv_sub$genes_collapsed)
            marg.lik0 <- dpois(genes_combo, n*cnv_sub$mut[match(names(genes_combo), cnv_sub$genes_collapsed)], log=TRUE)
            marg.lik1 <- dnbinom(genes_combo, cnv_sub$relative_risk[match(names(genes_combo), cnv_sub$genes_collapsed)]*beta.dn, 
              beta.dn/(beta.dn+n*cnv_sub$mut[match(names(genes_combo), cnv_sub$genes_collapsed)]), log=TRUE)
            lBF <- marg.lik1-marg.lik0
                lBF <- do.call(rbind, lapply(1:length(lBF), function(x) split_weight(lBF[x], info)))
        }
    }, cnv_use, n=n, del_dup_adj=del_dup_adj, prior=prior))

    sum_lBF <- by(raw_lBF$lBF, raw_lBF$gene, sum)
    BF <- rep(1, length(info$gene_gencodeV33))
        BF[match(names(sum_lBF), info$gene_gencodeV33)] <- exp(as.numeric(sum_lBF))
    return(BF)
}

################################################################################################
### Bayes factor calculation for case-control/inherited CNV contribution to a gene
## Input: 
# cnv_size_range: a vector of size range of CNV events to considered (# of constrained genes)
# cnv_use:  the set of cnvs to be fed into the model (for format refer to Supplementary Table 9 and run_TADA.R)
# n_case: the number of cases
# n_con: the number of controls (case-control analysis) or cses (inherited analysis)
# info: matrix that contains gene level information (Supplementary Table 8: gene names, constraint)
# prior: a vector of the prior risk for each gene
# del_dup_adj: the adjustment for duplications being in general less penetrant than deletions
# nu: a hyperparameter that controls nu0 and nu1 in determining prior distributions, passed to bayes.factor.cc
## Output:
# BF: a vector of Bayes factors of case-control/inherited evidence from CNVs for each gene for a given variant type
################################################################################################
BF_CC_CNV <- function(cnv_size_range, cnv_use, n_case, n_con, info, prior, del_dup_adj=1, nu=5000){

  n.cc.cnv <- data.frame(ca=n_case, cn=n_con)
  raw_lBF <- do.call(rbind, lapply(cnv_size_range, function(i, cnv_use){
    cnv_sub <- cnv_use[cnv_use$num_genes==i,]
    cnv_sub$relative_risk <- sapply(cnv_sub$genes, function(x) (sum(prior[match(x, info$gene_gencodeV33)])-1)/del_dup_adj+1)
    rho.ptv = nu * sum(cnv_sub$Affected_Status==1) / ((2*length(info$gene_gencodeV33)-i)*n.cc.cnv[2])

    gene_combos <- unique(cnv_sub$genes_collapsed)

    lBF_cc <- as.numeric(do.call(rbind, lapply(gene_combos, function(x){
      tot_lbf <- log(bayes.factor.cc(x.cc = data.frame(
        ca=sum(cnv_sub$genes_collapsed==x & cnv_sub$Affected_Status==2), 
        cn=sum(cnv_sub$genes_collapsed==x & cnv_sub$Affected_Status==1)), 
      n.cc = n.cc.cnv,
      gamma.cc=cnv_sub$relative_risk[cnv_sub$genes_collapsed==x][1], 
      rho1=as.numeric(rho.ptv), nu1=nu, 
      rho0=as.numeric(rho.ptv), nu0=nu))
      return(tot_lbf)
      })))
    names(lBF_cc) <- gene_combos
    lBF_cc <- do.call(rbind, lapply(1:length(lBF_cc), function(x) split_weight(lBF_cc[x], info)))
    }, cnv_use=cnv_use))
  sum_lBF <- by(raw_lBF$lBF, raw_lBF$gene, sum)
  BF <- rep(1, length(info$gene_gencodeV33))
      BF[match(names(sum_lBF), info$gene_gencodeV33)] <- exp(as.numeric(sum_lBF))
  return(BF)

}

################################################################################################
### Function to transform total BF evidence to FDR values
### Bayesian FDR control (PMID:19822692, Section2.3)
## Input: 
# BF: total Bayes factor matrix (each row is a gene, multiple columns allowed)
# pi0: estimated proportion of non-risk genes
## Output:
# FDR: a vector of FDR values corresponding to the rows of the input matrix
################################################################################################
Bayesian.FDR <- function(BF, pi0) {
  # order the BF in decreasing order, need to retain order to get results back in proper order 
  i.order=order(BF, decreasing = T)
  BF=BF[i.order]
  # convert BFs to PPA (posterior probability of alternative model)
  pi <- 1-pi0
  q <- pi*BF/(1-pi+pi*BF) # PPA
  q0 <- 1 - q # posterior probability of null model
  
  # the FDR at each PPA cutoff
  FDR=cumsum(q0)/(1:length(BF))
  
  # reorder to the original order
  FDR[i.order]=FDR
  
  return (FDR=FDR)
}

################################################################################################
### Helper function to label genes within a CNV with corresponding LOEUF, then subset by LOEUF threshold
## Input: 
# cnv:  a set of CNVs to be annotated (for format refer to Supplementary Table 9 and run_TADA.R)
# loeuf threshold: the threshold considered to be constrained and considered for modeling
# info: matrix that contains gene level information (Supplementary Table 8: gene names, constraint)
## Output:
# cnv: an annotated version of the input matrix, with fields for constrained genes (genes, genes_collapsed) and number of constrained genes (num_genes)
################################################################################################
processCNV <- function(cnv, loeuf_threshold=0.6, info){
  genes <- sapply(cnv$genes, function(x) unlist(str_split(x, ",")))
  genes_in <- lapply(genes, function(x) x[x %in% info$gene_gencodeV33])
  genes_loeuf <- lapply(genes_in, function(x) info$LOEUF[match(x, info$gene_gencodeV33)])

  genes_in <- lapply(1:length(genes), function(i) sort(genes_in[[i]][genes_loeuf[[i]]<=loeuf_threshold]))
  num_genes <- sapply(genes_in, length)
    names(num_genes) <- NULL
  cnv$genes <- genes_in
  cnv$num_genes <- num_genes
  cnv$genes_collapsed <- sapply(cnv$genes, function(x) paste0(x, collapse=","))
  return(cnv)
}

################################################################################################
### Helper function to distribute BF evidence inversely to LOEUF score for CNV evidence
## Input: 
# x: a vector of log Bayes factors evaluated for CNV events, with names labeled by a "," concatenated string of constituent genes
# info: matrix that contains gene level information (Supplementary Table 8: gene names, constraint)
## Output:
# evidence: a matrix containing the log Bayes (lBF) factor contributed to each gene(gene)from evaluated CNVs
################################################################################################
split_weight <- function(x, info){
    genes <- unlist(str_split(names(x), ","))
    oes <- info$LOEUF[match(genes, info$gene_gencodeV33)]
    weights <- sum(oes)/oes
        weights <- weights/sum(weights)
    return(evidence=data.frame(gene=genes, lBF=as.numeric(x)*weights))
}

################################################################################################
### Function to backtransform FDR to be on p-value scale
## Input: 
# x: a vector of FDRs
## Output:
# p: a vector of p-values
################################################################################################
q2p <- function(x){
    ec <- x
    ec <- sort(ec)
    ec <- ec*(1:length(ec))/(length(ec))/max(ec)    
    return(p=ec[names(x)])
}
