Configuring the pipeline
========================

Before you use the pipeline, you'll need to configure it. This configuration specifies how the scripts should run their commands and submit jobs to their cluster.

All the configuration details are in the config file ``train.cfg``.

The ``cluster`` option is used to determine how to parse the cluster's job status. For example, coyote uses ``qstat`` and the Broad uses ``bjobs``.

The ``bashrc`` option specifies a shell script to be sourced before any other commands are issued in a submitted job array. This allows us to specify that, say, a specific version of python be used. Canonical bashrc files are included in the ``bashrcs`` folder.