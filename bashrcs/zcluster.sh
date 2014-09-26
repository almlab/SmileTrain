# enable the module tool
#source /etc/profile.d/modules.sh

# load our preferred python version
#module add python/2.7.3
#imodule add biopython/1.5.9

export LD_LIBRARY_PATH=/usr/local/R/3.1.1/lib64/R/lib:/usr/local/openmpi/1.4.4/gcc444/lib:$LD_LIBRARY_PATH
export PATH=/home/mwalab/lancaste/bin:/usr/local/R/3.1.1/bin:$PATH
export PYTHONPATH=/home/mwalab/lancaste/lib:/home/mwalab/lancaste/lib/SmileTrain:/home/mwalab/lancaste/lib/dbOTUcaller:$PYTHONPATH

