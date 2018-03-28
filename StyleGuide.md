# Style Guide for SplineLib
In general the project sticks to the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html). This
guide only describes these project differs from Google's guide.

To check if the code is compliant with the Google C++ Style Guide [cpplint](https://github.com/cpplint/cpplint) can be
used. cpplint is a linter written in python. To install the linter see the guidelines on github.

## Formatting
### Line lenght
Google prescribes a [line length](https://google.github.io/styleguide/cppguide.html#Line_Length) of 80 characters. For
modern displays this is a very narrow restriction. Since variable and method names shall be descriptive, this limit is
reached very fast. Therefore, the line length is restricted to 120 characters. To account for this rule the change is
defined in the CPPLINT.cfg of this project.
