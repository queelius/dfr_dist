# Contributing to flexhaz

Thank you for your interest in contributing to flexhaz.

## Reporting Bugs

Open an issue at https://github.com/queelius/flexhaz/issues with:

- A minimal reproducible example
- Expected vs actual behavior
- Output of `sessionInfo()`

## Suggesting Features

Open an issue describing:

- The use case (what you're trying to do)
- How it relates to existing functionality
- Whether you'd be willing to implement it

## Submitting Changes

1. Fork the repository and create a branch from `main`
2. Make your changes following the conventions below
3. Add tests for new functionality (`tests/testthat/`)
4. Run `devtools::check()` — must pass with 0 errors and 0 warnings
5. Submit a pull request

## Development Conventions

- **Documentation**: roxygen2 comments on all exported functions. Run `devtools::document()` after changes. Never edit `NAMESPACE` or `man/` files by hand.
- **Testing**: testthat edition 3. Aim for high coverage. Run `devtools::test()`.
- **Style**: Standard R style. S3 methods named `generic.class`. Closures returned by distribution methods.
- **Dependencies**: Minimize new dependencies. Discuss in the issue before adding imports.

## Code of Conduct

This project follows the [Contributor Covenant](CODE_OF_CONDUCT.md). By participating, you agree to uphold this code.
