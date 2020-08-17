
## originally by Yulong Niu
## yulong.niu@hotmail.com

library('dcanr')

dcMethods()

data(sim102)
getConditionNames(sim102)

simdata <- getSimData(sim102, cond.name = 'UME6', full = FALSE)
emat <- simdata$emat
ume6_kd <- simdata$condition

z_scores <- dcScore(emat, ume6_kd, dc.method = 'zscore', cor.method = 'spearman')

raw_p <- dcTest(z_scores, emat, ume6_kd)
adj_p <- dcAdjust(raw_p, f = p.adjust, method = 'fdr')

library(igraph)

dcnet <- dcNetwork(z_scores, adj_p)
plot(dcnet, vertex.label = '')


adjmat <- as_adj(dcnet, sparse = FALSE)
edgedf <- as_data_frame(dcnet, what = 'edges')
print(head(edgedf))
