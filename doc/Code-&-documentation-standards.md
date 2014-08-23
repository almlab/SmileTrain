## Reading the in-code documentation
To look at docstrings in a particular script `script.py`, you can run `pydoc script`. Docstrings appear as red when the source is viewed on github.

## Coding standards
* Scripts and their functions have unit tests in the `test` folder. It is best that every function be tested; it's acceptable that just the few highest-level functions be tested. It's **essential** that there be *some* test, since this ensures all users of the pipeline that it hasn't been broken by any new changes.
    * Unit testing uses python's standard [unittest](https://docs.python.org/2/library/unittest.html).
* Indentation is 4 spaces. Tabs are not cool.
* Scripts that read input from the command line use argparse.

## Documentation standards
* Every function has a docstring.
    * Learn about [docstrings](http://legacy.python.org/dev/peps/pep-0257/).
    * Short, simple functions with one input and output can have one-line docstrings.
    * Functions that are more than about 5 lines, that have more than input, or that have an output whose type is not obvious should have a multiline docstring that specifies the purpose of the function, the expected types of inputs and outputs, the default values of any parameters, and the meanings of each parameter and output.