# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

recTree <- function(vars, obs, ntree, calculate_node_pred, X, Y, leftDaughter, rightDaughter, nodestatus, xbestsplit, nodepred, bestvar, inbag, varLevels, OOBtimes, localIncrements) {
    invisible(.Call('forestFloorStable_recTree', PACKAGE = 'forestFloorStable', vars, obs, ntree, calculate_node_pred, X, Y, leftDaughter, rightDaughter, nodestatus, xbestsplit, nodepred, bestvar, inbag, varLevels, OOBtimes, localIncrements))
}

