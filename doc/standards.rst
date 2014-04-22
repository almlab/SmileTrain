Coding & documentation standards
=======================================

Reading the documentation
-------------------------

Most of the documentation that is in the source code is automatically included in these documents. To look at docstrings in a particular script ``script.py``, you can run ``pydoc script``.

Coding standards
----------------
* Scripts and their functions have unit tests in ``tests.py``. It is best that every function be tested; it's acceptable that just the few highest-level functions be tested. It's **essential** that there be *some* test, since this ensures all users of the pipeline that it hasn't been broken by any new changes.
* Indentation is 4 spaces. Tabs are not cool.
* Scripts that read input from the command line use argparse.

Documentation standards
-----------------------
* Every script has a docstring. This should provide a description of the script's purpose and a quick overview of what its inputs and outputs should look like. For example:
* Every function has a docstring. Short, simple functions with one input and output can have one-line docstrings. Functions that are more than about 5 lines, that have more than input, or that have an output whose type is not obvious should have a multiline docstring that specifies the purpose of the function, the expected types of inputs and outputs, the default values of any parameters, and the meanings of each parameter and output. The format for multiline function docstrings is::

    '''
    Description of the function, along with its intended use; say, releasing passenger pigeons.
    
    destinations : list of strings
        destination cities for the birds
    direction : string ('north', 'south', 'east', 'west'; default 'north')
        cardinal direction that the bird will fly in
    n_birds : int (default 1)
        number of birds released
        
    returns : int
        number of birds that made it back alive
    '''
