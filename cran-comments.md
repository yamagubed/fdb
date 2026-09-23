## Release preparation status

This is a local preparation note, not a claim of completed CRAN checks.
No submission has been made.

The current maintainer email and repository URLs were retained at the
maintainer's request. The example.com email is a submission blocker until
a functioning maintainer address is supplied. Verify the public URLs and
package-name availability before uploading.

Source-level regression tests and documentation checks were run in an
isolated webR environment. Native Windows/macOS/Linux R CMD check,
native PSOCK integration tests, and R-devel checks must be completed before
replacing this note with an actual check-results summary.

Large design simulations are excluded from automatic vignette execution.
Parallel execution defaults to two workers. Native parallel integration
tests are enabled by FDB_RUN_PARALLEL_TESTS=true and use two workers.
