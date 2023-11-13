The following is information than may be useful for contributors. The information is advisory in nature.

All communication is held in [issues](https://github.com/stencillogic/num-bigfloat/issues) and [PRs](https://github.com/stencillogic/num-bigfloat/pulls).

There are several ways to contribute to the project. The list includes, but is not limited by:

 - opening issues, participating in discussions
 - implementing features, or fixing bugs which are defined by current open issues, opening pull requests
 - writing documentation, or articles about the library
 - developing and improving algorithms used in the library

Rust is the primary programming language of the library. The following are reccomendations to follow with regard to the coding standards:

 - use of unsafe code must be kept at a minimum; any use of unsafe must be justified.
 - use of `unwrap()`, `expect()`, any unhandled errors are only acceptable when it was proved that they will never appear (i.e. there are checks that guarantee the code will never panic or return error).
 - new code must be covered by regression tests.
 - code must be formatted with [rustfmt](https://rust-lang.github.io/rustfmt/).
 - all public items (functions, structs, etc.) must be documented.
 - any pull request must pass all tests before it can be merged.
 - public api should follow [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/).
 - library's documentation should be kept up to date.
 - common sense and proper reasoning supersede any written rule.
