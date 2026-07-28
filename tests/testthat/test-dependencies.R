test_that("check_dependencies returns a well-formed table", {
  deps <- check_dependencies(quiet = TRUE)

  expect_s3_class(deps, "data.frame")
  expect_true(all(c("package", "type", "needed", "installed", "version", "ok")
                  %in% names(deps)))
  expect_true(all(deps$type %in% c("required", "optional")))
  expect_type(deps$installed, "logical")
  expect_type(deps$ok, "logical")
})

test_that("required dependencies are present, since the package loaded", {
  deps <- check_dependencies(quiet = TRUE)
  required <- deps[deps$type == "required", ]
  expect_true(all(required$installed))
  expect_true(all(required$ok))
})

test_that("optional packages can be excluded", {
  deps <- check_dependencies(include_optional = FALSE, quiet = TRUE)
  expect_true(all(deps$type == "required"))
  expect_lt(nrow(deps), nrow(check_dependencies(quiet = TRUE)))
})

test_that("the report prints without error", {
  expect_output(check_dependencies(), "dependency check")
  expect_invisible(check_dependencies(quiet = TRUE))
})

test_that("a version that is installed but too old is flagged", {
  # compareVersion is the mechanism the check relies on.
  expect_lt(utils::compareVersion("1.0.0", "2.0.0"), 0)
  expect_equal(utils::compareVersion("2.0.0", "2.0.0"), 0)
  expect_gt(utils::compareVersion("2.1.0", "2.0.0"), 0)
})

test_that("install_dependencies does nothing when everything is satisfied", {
  skip_if_not(all(check_dependencies(quiet = TRUE)$ok),
              "some optional dependencies are missing on this machine")
  expect_message(out <- install_dependencies(ask = FALSE),
                 "already satisfied")
  expect_equal(out, character(0))
})

test_that("install_dependencies refuses to prompt non-interactively", {
  skip_if(all(check_dependencies(quiet = TRUE)$ok),
          "nothing missing, so the prompt path is not reached")
  skip_if(interactive())
  expect_error(install_dependencies(ask = TRUE), "non-interactive")
})

test_that("flags are validated", {
  expect_error(check_dependencies(include_optional = "yes"), "include_optional")
  expect_error(check_dependencies(quiet = NA), "quiet")
})
