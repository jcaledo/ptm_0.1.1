library(ptm)
context("PTM Properties Representations")

## ---------------------------------------------- ##
#             Testing plot.ptm                     #
## ---------------------------------------------- ##

test_that("ptm.plot() works properly",{

  skip_on_cran()
  skip_on_travis()

  a <- ptm.plot("Q01196", property = 'sasa', ptm = 'p',
                dssp = 'compute', window = 1, sdata = TRUE) # ---------- single chain
  b <- ptm.plot("Q01196", pdb = "1H9D.C", property = 'acc', ptm = 'p',
                dssp = 'compute', window = 1, sdata = TRUE) # ------ non-protein chains
  c <- suppressWarnings(ptm.plot('G3SB67', property = 'eiip', window = 10, ptm = 'all')) # --- non-ptm sites
  d <- ptm.plot('P23246', property = 'acc', ptm = 'all') # ---  PDB only 42 % coverage

  expect_is(a, 'character')
  expect_true(grepl("Work done.", a))

  expect_is(b, 'character')
  expect_true(grepl("Work done.", b))

  expect_is(c, 'character')
  expect_true(grepl("Work done.", c))
})


# test_that("plot.ptm() works properly",{
#
#   skip_on_cran()
#   skip_on_travis()
#
#   a <- plot.ptm("P09803", property = 'eiip', window = 10, ptm = c('meto', 'p', 'ac'))
#   b <- plot.ptm("P09803", property = 'dpx', ptm = c('meto', 'p', 'ac'))
#   c <- plot.ptm("P09803", pdb = '1I7W.D', property = 'acc', ptm = c('meto', 'p', 'ac'))
#
#   expect_is(a, 'character') # --------------- length uniprot seq longer than the pdb seq
#   expect_true(grepl("Work done.", a))
#
#   expect_is(b, 'character') # -------------- length pdb shorter than seq uniprot
#   expect_true(grepl("Work done.", b))
#
#   expect_is(c, 'character') # -------------- length pdb with gaps
#   expect_true(grepl("Work done.", c))
# })

# test_that("plot.ptm() works properly",{
#
#   skip_on_cran()
#   skip_on_travis()
#
#   a <- plot.ptm("P18031", property = 'dpx', window = 1, ptm = 'all') # --- PTMs mostly on surface
#   b <- plot.ptm("P84022", property = 'dpx', window = 1, ptm = 'all') # -- PTM only on surface
#   # c <- plot.ptm("P02511", property = 'dpx', window = 1, dssp = 'mkdssp', ptm = 'all') # --- Disabling
#
#   expect_is(a, 'character')
#   expect_true(grepl("Work done.", a))
#
#   expect_is(b, 'character')
#   expect_true(grepl("Work done.", b))
# })


test_that("find.aaindex() works properly",{

  a <- find.aaindex('mutability')
  b <- find.aaindex('Kyte-Doolittle')
  c <- find.aaindex('argos')
  d <- find.aaindex('tontuna')

  expect_is(a, 'integer')
  expect_is(names(a), 'character')
  expect_equivalent(a, c(65, 135))
  expect_equal(names(a), c("DAYM780201", "JOND920102"))

  expect_is(b, 'integer')
  expect_is(names(b), 'character')
  expect_equivalent(b, c(151, 494))
  expect_equal(names(b), c("KYTJ820101", "JURD980101"))

  expect_is(c, 'integer')
  expect_is(names(c), 'character')
  expect_equivalent(c, c(2, 3, 4))
  expect_equal(names(c), c("ARGP820101", "ARGP820102", "ARGP820103"))

  expect_is(d, 'integer')
  expect_equal(length(d), 0)
})
