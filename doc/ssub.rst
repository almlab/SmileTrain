ssub
====

`ssub` is a tool for submitting many jobs. Some clusters have restrictions on how jobs may be submitted, and ssub deals with those restrictions.

What ssub does
--------------

Clusters use a queuing system. Coyote uses PBS. A queueing system allows multiple users to submit multiple command-line instructions or *tasks* to the cluster. In the simplest case, each tasks is a single *job* that takes up a single CPU. However, a user might want to submit many tasks while still allowing other users to use the cluster. One job might have many tasks in serial, but this can be slow, since each task proceeds in sequence.

A user can instead submit a *job array*, which is a list of jobs that get split up into a set number of CPUs. If there are 20 tasks and a user wants to use 2 CPUs, the job array will assign 10 tasks to one CPU and 10 to the other.

Certain cluster have restrictions, however. Coyote only allows user to submit job arrays with 500 jobs, regardless of the number CPUs required. ssub creates no more than 500 jobs, each composed of many tasks, and joins them in a job array. In this way, the cluster's requirements are met and the user has the fewest number of computational tasks per job.

Classes & Functions
-------------------

.. automodule:: ssub
	:members:
	:undoc-members: