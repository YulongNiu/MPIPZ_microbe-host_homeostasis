
if (!require(vegan)) install.packages('vegan')
library(vegan)

variability_table <- function(cca){

        chi <- c(cca$tot.chi,
                       cca$CCA$tot.chi, cca$CA$tot.chi)
        variability_table <- cbind(chi, chi/chi[1])
        colnames(variability_table) <- c("inertia", "proportion")
        rownames(variability_table) <- c("total", "constrained", "unconstrained")
        return(variability_table)

}

cap_var_props <- function(cca){

        eig_tot <- sum(cca$CCA$eig)
        var_propdf <- cca$CCA$eig/eig_tot
        return(var_propdf)
}

pca_var_props <- function(cca){

        eig_tot <- sum(cca$CA$eig)
        var_propdf <- cca$CA$eig/eig_tot
        return(var_propdf)
}

cca_ci <- function(cca, permutations=5000){

        var_tbl <- variability_table(cca)
        p <- permutest(cca, permutations=permutations)
        ci <- quantile(p$F.perm, c(.05,.95))*p$chi[1]/var_tbl["total", "inertia"]
        return(ci)

}
