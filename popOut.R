library(data.table)

assign.to.quants <- function( y, n ){
    q <- quantile( y, 1:(n-1)/n, na.rm=TRUE )
    y.cat <- sapply( lapply( y, '>=', q ), sum )
    return(y.cat)
}

prs.test <- function( prs, y, K=0.01 ){# Test for residuals mean=0 in the tail(s)
    
    #print(prs) 
    
    fit <- lm( prs ~ y )
    # Test upper tail "regression to mean"
    T <- quantile( y, 1-K, na.rm=TRUE )
    ptr <- which(y>T)
    test <- t.test( fit$residual[ptr] )#, alternative='less' )
    p2 <- test$p.value
    effect2 <- test$estimate
    ci2 <- test$conf.int
    # Test lower tail "regression to mean"
    T <- quantile( y, K, na.rm=TRUE )
    ptr <- which(y<T)
    test <- t.test( fit$residual[ptr] )#, alternative='greater' )
    p1 <- test$p.value
    effect1 <- test$estimate
    ci1 <- test$conf.int
    ms = sd(prs) 
    ret <- c( effect1, ci1, p1, effect2, ci2, p2, ms, effect1/ms, ci1/ms,effect2/ms, ci2/ms)
    names(ret) <- c( 'effect.lower', 'ci.lower.lower', 'ci.lower.upper', 'p.lower',
                    'effect.upper',' ci.upper.lower',' ci.upper.upper', 'p.upper' ,'sd',
                    're.lower','ri.lower.lower','ri.lower.upper','re.upper','ri.upper.lower','ri.upper.upper')            
    
    return(ret)
}


prs.test.qc <- function( prs, y, n.q=100, K=0.1, K1=0.9 ){
    q <- assign.to.quants( y, n.q )
    
    
    
    
    ptr <- which( n.q*K <= q&q < n.q*K1 )
    fit <- lm( prs ~ y, subset=ptr )
    

    
    p <- vector()
    
    for( j in unique(q[ptr]) ){
        ptr1 <- which( q[ptr]==j )
        p <- c( p, t.test( fit$residuals[ptr1] )$p.value )
    }
    qc.test <- ifelse( min(p) < 0.05/length(unique(q[ptr])), 'FAIL', 'PASS' )
        
    return(c( qc.test, min(p)*length(unique(q[ptr])) ))
}





prs.test.emp <- function( prs, y, n.q=100 ){
    stat <- vector()
    fit <- lm( prs ~ y )


    q <- assign.to.quants( y, n.q )
    
    
    for( i in 1:n.q ){
        ptr <- which(q==(i-1))
        test <- t.test( fit$residual[ptr] )
        stat[i] <- test$statistic
        if (i < 2){ 
            print(i)
            print(stat[i])  
        }
    }
    r1 <- rank(-stat)
    r2 <- rank(stat)
   
    p <- c( r1[1]/n.q, r2[100]/n.q )
    
    return(p)
}

args <- commandArgs(trailingOnly=TRUE)
prs.file <- args[1]
pheno.file <- args[2]

#prs <- fread(prs.file)
prs <- na.omit(fread(prs.file))


#pheno <- fread(pheno.file)
pheno <- na.omit(fread(pheno.file))






ids <- intersect( prs$IID, pheno$IID )


size <- length(ids) 




prs <- prs[match( ids, prs$IID ),]
pheno <- pheno[match( ids, pheno$IID ),]





test <- prs.test( prs$PRS, pheno$Phenotype )

test.empirical <- prs.test.emp( prs$PRS, pheno$Phenotype )

qc <- prs.test.qc( prs$PRS, pheno$Phenotype )

t_name = strsplit(pheno.file, split=".target") 

out <- data.frame( t_name, size, t(c(test, test.empirical, as.numeric(qc[2]))), qc[1] )


#print(out) 


names(out) <- c( 'file', 'size', 
                'effect.lower', 'ci.lower.lower', 'ci.lower.upper', 'p.lower',
                'effect.upper','ci.upper.lower','ci.upper.upper', 'p.upper',
                'sd','re.lower','ri.lower.lower','ri.lower.upper','re.upper','ri.upper.lower','ri.upper.upper',       
                'p.emp.lower','p.emp.upper','p.qc','QC')
    



out <- as.data.frame(t(out))



print( out )

#pheno.file='file'
#prs=data.frame(1:1000,rnorm(1000))
#pheno=data.frame(1:1000,rnorm(1000))
#colnames(pheno)=c('IID','Phenotype')
#colnames(prs)=c('IID','PRS')
#ids=1:1000















