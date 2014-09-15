SmileTrain
==========

## Python requirements
SmileTrain is developed against Python 2.7.3.

## Documentation
Documents in the github wiki (in the right-hand panel). A (possibly old) copy is in the `doc/` folder.

Scripts are annotated with docstrings that can be read using pydoc. To view documentation for script.py, run `pydoc script`.

## Tests
Unit tests are in the `test` folder. You can run the tests from the top directory by running `python -m unittest discover`.

## Usage
* Enter your user information in the config file before use.
* Make sure your .bash\_rc loads python 2.7. Otherwise, argparse will die. At the start of my .bash\_rc, I have
    - source /etc/profile.d/modules.sh
    - module add python/2.7.3
* Run the tests.
* Run otu\_caller.py with the command line option --dry\_run before you actually submit any jobs. Check that the output commands look like what you actually want to submit to the cluster!
