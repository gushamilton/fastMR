# The "worthwhile threads" rule keeps tiny inputs on one worker. The existing
# thread-equivalence tests use tiny inputs on purpose, so force the parallel
# code paths for the whole run; test-thread-cap.R restores the real rule.
previous_work_scale <- fastMR:::fastmr_set_work_scale_native(0)
withr::defer(fastMR:::fastmr_set_work_scale_native(previous_work_scale),
             testthat::teardown_env())
