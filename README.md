smile_train
===========

# Documentation
The full documentation is in the `doc` subfolder. To compile on Linux or Mac,

    make html

On Windows, you need to use make.bat somehow.

Then you can open `index.html` and read away!

Scripts are annotated with docstrings that can be read using pydoc. To view documentation for script.py, run
    pydoc script

# Dependencies
* Python 2.7, numpy, pandas

# Usage
* Enter your user information in the config file before use.
* Make sure your .bash_rc loads python 2.7. Otherwise, argparse will die. At the start of my .bash_rc, I have
    - source /etc/profile.d/modules.sh
    - module add python/2.7.3
* Run the tests.
* Run otu_caller.py with the command line option --dry_run before you actually submit any jobs. Check that the output commands look like what you actually want to submit to the cluster!

# Code standards
* Indentation is 4 spaces
* Functions have a docstring. Functions with multiple inputs or multiple outputs should have a multiline docstring specifying the type of each input and output value.
* Scripts have a docstring.
* Scripts that read input from the command line use argparse.
