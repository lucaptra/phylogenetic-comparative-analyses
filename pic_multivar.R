dyn.load(paste("pic_multivar", .Platform$dynlib.ext, sep = ""))

pic_multivar <- function (x, phy, scaled = TRUE, var.contrasts = FALSE, rescaled.tree = FALSE)
{
    if (!inherits(phy, "phylo")) {
        stop("object 'phy' is not of class \"phylo\"")}
    
    if (is.null(phy$edge.length)) {
        stop("your tree has no branch lengths")}
    
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    
    if (nb.node != nb.tip - 1) {
        stop("'phy' is not rooted and fully dichotomous")}
    
    if (any(is.na(x))) {
        stop("missing data in 'x': you may consider removing the species with missing data from your tree with the function 'drop.tip'.")}
    
    phy <- reorder(phy, "postorder")
    phenotype <- matrix(0, length(numeric(nb.tip + nb.node)), length(x[1, ]))
    
    if (is.null(row.names(x))) {
        phenotype[1:nb.tip, ] <- x
    } else {
        if (all(row.names(x) %in% phy$tip.label)) {
            phenotype[1:nb.tip, ] <- matrix(sapply(x[row.names(x) == phy$tip.label, ], as.numeric), length(x[row.names(x) == phy$tip.label, 1]), length(x[row.names(x) == phy$tip.label[1], ]))
        } else {
            phenotype[1:nb.tip, ] <- x
            warning("the names of argument 'x' and the tip labels of the tree did not match: the former were ignored in the analysis.")}
        }
    
    ans <- .C("C_pic_multivar", as.integer(nb.tip), as.integer(nb.node), as.integer(phy$edge[, 1]), as.integer(phy$edge[, 2]), as.double(phy$edge.length), as.double(phenotype), as.integer(dim(phenotype)[1]), as.integer(dim(phenotype)[2]), double(nb.node), double(nb.node), as.integer(var.contrasts), as.integer(scaled))
    
    contr <- ans[[9]]
    
    lbls <- if (is.null(phy$node.label)) {
        as.character(1:nb.node + nb.tip)
       } else {phy$node.label}
    
    if (var.contrasts) {
        contr <- cbind(contr, ans[[10]])
        dimnames(contr) <- list(lbls, c("contrasts", "variance"))
    } else {names(contr) <- lbls}
    
    if (rescaled.tree) {
        phy$edge.length <- ans[[5]]
        contr <- list(contr = contr, rescaled.tree = phy)
    }
    
    contr
}