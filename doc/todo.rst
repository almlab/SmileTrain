Coding & documentation to do
============================

chimeras
--------
* Nothing here has been refactored yet. It might not even work.
* Shera should be implemented here.

otu_caller.py
-------------
* Put a format check before doing quality filtering. If the format is wrong, then usearch dies and ruins the output.
* Implement a more effective dry run. Right now the dry run checks for the existence of some files and throws errors. It would be nicer if it just told you where it was goign to check for file existence and not stop the dry run.

intersect
---------
* Does fastq_ids belong here? Maybe should go to utils if someone else uses it

ref_gg
------
* Something is screwy with uclust's measurement of % identity. I've seen alignments like 250I250M250I that get 99.6% identity and alignments like 250I200MD49M250I that get 100%.

dbotus
------
* Integrate SPP's dbOTU calling code.

ssub
----
* Refactor and remove old code. Decide if this ssub should be entirely a module for the otu caller. Otherwise, it needs to be able to be called as a script. Otherwise, ssub and the otu caller will have to be separate projects that depend on each other, which sounds a little dangerous.
* Put in a check for available clusters.
* Add functionality & tests for broad.
    - qstat equivalent
    - bashrc
* Double-check the way job names are stored. Should keep integers and the [], since that's what you can use with ``qdel``.

doc
---
* Collect some more SPP wisdom.

Diagnostics
-----------
* Figure out where to put diagnostics. A separate folder?
* Cook up some diagnostics:
    - Are the barcodes pointing the right way? Not rev comp'd?
    - How many barcodes are good? How many have one error?
    - What fraction of reads got thrown out at each step? Scripts give some diagnostic error message after finishing?
    - Keep a log?
    
Strategic
---------
* Move away from a nohup/screen kind of usage.
    - Log files would help a lot with that!
* Put more things in configuration files?
    - e.g., all the command line options
    - this would serve as a record for how the run went
* Keep logs of results?
    - in the log, keep track of the hash of the output file to see if you have the same one or if it's changed
* Refactor the message and warning functions.
    - Errors should raise exceptions
    - Warnings should raise warnings?
    - Message should interact with some log?
* How do we keep track of errors?
    - Should be some bug tracker
    - Triaging