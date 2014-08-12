
estimate_ssGSEA_BB<-function(exprs=NULL,SI_geneset=NULL){
    m <- data.matrix(exprs)
    gene_names <- rownames(exprs)
    sample_names <- colnames(exprs)
    sample_count <- length(m[1, ])
    gene_count <- length(m[, 1])
    for (j in 1:sample_count) {
        m[, j] <- rank(m[, j], ties.method = "average")
    }
    m <- 10000 * m/gene_count
    gs <- SI_geneset
    N.gs <- length(SI_geneset)
    gs.names <- names(SI_geneset)
    score.matrix <- matrix(0, nrow = N.gs, ncol = sample_count)
    for (gs.i in 1:N.gs) {
        gene.set <- gs[[gs.i]]
        gene.overlap <- intersect(gene.set, gene_names)
        print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", 
            length(gene.overlap)))
        if (length(gene.overlap) == 0) {
            score.matrix[gs.i, ] <- rep(NA, sample_count)
            next
        }
        else {
            ES.vector <- vector(length = sample_count)
            for (S.index in 1:sample_count) {
                gene.list <- order(m[, S.index], decreasing = TRUE)
                gene.set2 <- match(gene.overlap, gene_names)
                correl.vector <- m[gene.list, S.index]
                TAG <- sign(match(gene.list, gene.set2, nomatch = 0))
                no.TAG <- 1 - TAG
                N <- length(gene.list)
                Nh <- length(gene.set2)
                Nm <- N - Nh
                correl.vector <- abs(correl.vector)^0.25
                sum.correl <- sum(correl.vector[TAG == 1])
                P0 <- no.TAG/Nm
                F0 <- cumsum(P0)
                Pn <- TAG * correl.vector/sum.correl
                Fn <- cumsum(Pn)
                RES <- Fn - F0
                max.ES <- max(RES)
                min.ES <- min(RES)
                if (max.ES > -min.ES) {
                  arg.ES <- which.max(RES)
                } else {
                  arg.ES <- which.min(RES)
                }
                ES <- sum(RES)
                EnrichmentScore <- list(ES = ES, arg.ES = arg.ES, 
                  RES = RES, indicator = TAG)
                ES.vector[S.index] <- EnrichmentScore$ES
            }
            score.matrix[gs.i, ] <- ES.vector
        }
    }
    score.data <- data.frame(score.matrix)
    colnames(score.data) <- sample_names
    rownames(score.data) <- gs.names
    return(score.data)
}