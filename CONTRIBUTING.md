# Contributing to parcutils

This project uses standard R package tooling. To run local package checks or to build the package, you must have R installed and available on your `PATH`.

The typical workflow is:

```sh
R CMD build .
R CMD check parcutils_*tar.gz
```

Our continuous integration setup installs R using [r-lib/actions/setup-r](https://github.com/r-lib/actions/tree/master/setup-r). Ensure you have a compatible R installation when running these commands locally to avoid check failures.
