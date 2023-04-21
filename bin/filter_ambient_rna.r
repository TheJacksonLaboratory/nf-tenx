#!/usr/bin/env Rscript

library(SoupX)
library(DropletUtils)

sc = load10X('.')
sc = autoEstCont(sc)
out = adjustCounts(sc)

write10XCounts('./soupx_filt_mat')