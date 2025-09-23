get_nodes <- function(adj, num) {
    N <- length(num)
    nn <- num
    N_edges <- length(adj) / 2
    node1 <- vector(mode = "numeric", length = N_edges)
    node2 <- vector(mode = "numeric", length = N_edges)
    iAdj <- 0
    iEdge <- 0
    for (i in 1:N) {
        for (j in 1:nn[i]) {
            iAdj <- iAdj + 1
            if (i < adj[iAdj]) {
                iEdge <- iEdge + 1
                node1[iEdge] <- i
                node2[iEdge] <- adj[iAdj]
            }
        }
    }
    return(list("N" = N, "N_edges" = N_edges,
                "node1" = node1, "node2" = node2))
}
