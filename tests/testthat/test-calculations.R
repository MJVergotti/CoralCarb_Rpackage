test_that("Temperature conversion functions correctly", {
  expect_equal(TCel_to_Kelv(0), 273.15)
  expect_equal(TCel_to_Kelv(25), 298.15)
})

test_that("Salinity and constant equations return valid positive values", {
  expect_true(Btot(35) > 0)
  expect_true(Ksp_arag(25, 35) > 0)
  expect_true(Kb_sw(25, 35) > 0)
})

test_that("pHsw_coral matches expected species calibration values", {
  # Cladocora caespitosa calibration test
  # Delta_pH = 4.7974 - 0.521 * pHsw
  # pHcf - pHsw = 4.7974 - 0.521 * pHsw => pHcf - 4.7974 = 0.479 * pHsw
  # For pHcf = 8.7, pHsw should equal ~8.14739
  expect_equal(pHsw_coral(8.7, species = "Cladocora caespitosa"), 8.14739, tolerance = 1e-4)
})

test_that("coralCF returns a complete data.frame with expected structure and values", {
  res <- coralCF(TCel = 25, Sal = 35, d11B = 24.43, BCa = 664.8, species = "Cladocora caespitosa")
  
  # Structural checks
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  
  # Confirm required columns exist
  required_cols <- c("TCel", "Sal", "pHcf", "pHsw", "borate_cf", "carbonate_cf", "DICcf", "arg_sat_cf")
  expect_true(all(required_cols %in% colnames(res)))
  
  # Numerical sanity checks
  expect_gt(res$pHcf, 8.0)
  expect_gt(res$arg_sat_cf, 1.0)
})

test_that("coralCF handles vector inputs correctly", {
  res_vec <- coralCF(
    TCel = c(25, 26), 
    Sal = c(35, 35), 
    d11B = c(24.43, 24.80), 
    BCa = c(664.8, 670.0)
  )
  
  expect_equal(nrow(res_vec), 2)
  expect_equal(res_vec$TCel, c(25, 26))
})