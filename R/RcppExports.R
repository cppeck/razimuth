# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

priorUpts <- function(gridpts, obspt, radius) {
    .Call('_razimuth_priorUpts', PACKAGE = 'razimuth', gridpts, obspt, radius)
}

vhf_cc_mu <- function(cp, thetalist, thetatmp, distid) {
    .Call('_razimuth_vhf_cc_mu', PACKAGE = 'razimuth', cp, thetalist, thetatmp, distid)
}

pmode_cpp <- function(x, ux) {
    .Call('_razimuth_pmode_cpp', PACKAGE = 'razimuth', x, ux)
}

