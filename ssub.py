#!/usr/bin/env python

import argparse, os, re, select, stat, subprocess, sys, tempfile, time, re, time, itertools
import xml.etree.ElementTree as ET
from util import *

'''
ssub is a simple job submission script

to use it as a python library:
  import ssub
  ssub.args.n = 100
  ssub.args.q = "hour"
  commands = get_commands()
  job_ids = submit_array(ssub.args, commands)
  wait_for_jobs(job_ids)
  print "done"

to create and run pipelines:
  import ssub
  ssub.args.q = "week"
  ssub.args.G = "broadfolk"
  commands1 = run_blast()
  commands2 = parse_blast()
  pipeline = [commands1, commands2]
  submit_pipeline(args, pipeline)
  print "done"

'''

def coyote_parse(qstat_output, my_username):
    '''
    Parse output from qstat -x, looking for the jobs that have the user as owner

    qstat_output : string
        output from qstat -x
    my_username : string
        kerberos username

    Returns : list of job IDs, each a string of integers with an optional []
    '''

    root = ET.fromstring(qstat_output)
    job_ids = []
    for child in root:
        # raw username contains whoever@wiley.coyote.etc.etc
        # get just the part that precedes the @
        raw_username = child.find('Job_Owner').text
        username = re.match('(.+)@', raw_username).group(1)

        if username == my_username:
            # job id contains 12345.wiley.coyote.etc.etc
            # we want just the integers at the beginning
            raw_job_id = child.find('Job_Id').text
            job_id = re.match('\d+(\[\])?', raw_job_id).group()
            job_ids.append(job_id)

    return job_ids

def zcluster_parse(qstat_output, my_username):
    '''parse output from `qstat -xml`'''
    root = ET.fromstring(qstat_output)
    job_ids = []

    # root should have two children, either queue_info or job_info
    for qinfo_or_jinfo in root:
        # those children should have elements job_list
        for job_list in qinfo_or_jinfo:
            # job_list node should contain owner and job number info
            username = job_list.find('JB_owner').text

            if username == my_username:
                job_id = job_list.find('JB_job_number').text
                job_ids.append(job_id)

    return job_ids


class Ssub():
    def __init__(self, username, cluster, queue, tmp_dir, bashrc, n_cpus=1, header='#!/bin/bash', l=200, m=4, group='', io=''):
        # initialize cluster parameters
        self.cluster = cluster
        self.username = username
        self.tmp_dir = tmp_dir
        self.q = queue

        self.header = header + "\n"
        
        self.l = l  # slot limit?
        self.n_cpus = n_cpus
        self.m = m  # gb of memeory
        
        self.G = group  # broad specific?
        self.io = io    # broad specific?

        self.bashrc = bashrc
        self.source_line = 'source %s\n' % self.bashrc
        
        # broad parameters
        if self.cluster == 'broad':
            raise RuntimeError("broad may not be supported")
            # swo> stat_cmd should be a list of words passed to terminal
            self.submit_cmd = 'bsub'
            self.stat_cmd = 'bjobs -w'
            self.parse_job = lambda x: re.search('Job <(\d+)>', x).group(1)
            # swo> should parse the status!
            self.parse_stats = lambda x: 0
        
        # coyote parameters
        elif self.cluster == 'coyote':
            self.submit_cmd = 'qsub'
            self.stat_cmd = ['qstat', '-x']
            # parse: match a string of integers then either [] or nothing
            self.parse_job = lambda x: re.match('\d+(\[\])?', x).group()
            self.parse_status = lambda x: coyote_parse(x, username)

        elif self.cluster == 'zcluster':
            self.submit_cmd = 'qsub' 
            self.stat_cmd = ['qstat', '-xml']
            self.parse_job = lambda x: re.match('(\d+)\.', x.split()[2]).group(1)
	    self.parse_status = lambda x: zcluster_parse(x, username)
        
        else:
            raise RuntimeError('unrecognized cluster %s' %(cluster))
    
    def __repr__(self):
        return '\ncluster: %s\nusername: %s\ntmp_dir: %s\nn: %d\nm: %d\nq: %s\nG: %s\nio: %s\n' %(self.cluster, self.username, self.tmp_dir, self.n_cpus, self.m, self.q, self.G, self.io)
    
    def __str__(self):
        a = self.__repr__()
        a = re.sub('\n', '\n  ', a)
        return a
    
    def mktemp(self, suffix='.tmp'):
        '''make temporary file and return [fh, fn]'''
        fh, fn = tempfile.mkstemp(suffix='.sh', dir=self.tmp_dir)
        os.close(fh)
        fh = open(fn, 'w')
        fh.write(self.header)
        return fh, fn
    
    def job_status(self):
        '''call the command specified by the words in stat_cmd'''
        return subprocess.check_output(self.stat_cmd)
    
    def jobs_finished(self, job_ids):
        '''
        Check if jobs are finished

        job_ids : list of strings
            ids to search for in the output of job status command

        Returns : bool
            were any of the ids found in the job status output?
        '''

        status = self.job_status()
        running_jobs = self.parse_status(status)

        for job in running_jobs:
            if job in job_ids:
                return False
            
        message('jobs completed\ttime: %s' % time.strftime("%d %b %H:%M", time.localtime()), indent=6)
        return True
    
    def write_jobs(self, commands):
        '''write job scripts from a list of commands'''
        
        # initialize output files
        fhs, fns = zip(*[self.mktemp(suffix='.sh') for i in range(min(self.n_cpus, len(commands)))])
        fhs_cycle = itertools.cycle(fhs)
        
        # write commands to file
        for command in commands:
            fh = fhs_cycle.next()
            fh.write('%s\n' %(command))
        
        # close all filehandles
        for fh in fhs:
            fh.close()
        
        # make executable and print message
        for fn in fns:
            os.chmod(fn, stat.S_IRWXU)
            message('Writing job %s' %(fn), indent=4)
        
        return fns
    
    def submit_jobs(self, fns):
        '''submit jobs to the cluster'''
        job_ids = []
        for fn in fns:
            if self.cluster == 'broad':
                process = subprocess.Popen(['%s < %s' %(self.submit_cmd, fn)], stdout = subprocess.PIPE, shell=True)
            elif self.cluster in ['coyote', 'zcluster']:
                process = subprocess.Popen(['%s %s' %(self.submit_cmd, fn)], stdout = subprocess.PIPE, shell=True)
            else:
                raise RuntimeError("trying to submitting job on unsupported cluster")

            [out, error] = process.communicate()

            try:
                job_id = self.parse_job(out)
            except:
                raise RuntimeError("could not parse output='%s' error='%s' with job parser" %(out, error))

            job_ids.append(job_id)
            message('Submitting job %s' %(fn), indent=4)
            message('job ID: %s\ttime: %s' %(job_id, time.strftime("%d %b %H:%M", time.localtime())), indent=6)
        return job_ids
    
    def write_LSF_array(self, fns):
        '''write an LSF job array from args and filenames'''
        
        # initialize output file
        fh, array_fn = self.mktemp(suffix='.sh')
        array_fn = os.path.abspath(array_fn)
        
        # write header
        fh.write('#BSUB -J "job[1-%d]%%%s"\n' %(len(fns), self.l))
        fh.write('#BSUB -e %s.e.%%I\n' %(array_fn))
        fh.write('#BSUB -o %s.o.%%I\n' %(array_fn))
        fh.write('#BSUB -q %s\n' %(self.q))
        fh.write('#BSUB -G %s\n' %(self.G))
        fh.write('#BSUB -R "rusage[mem=%s:argon_io=%s]"\n' %(self.m, self.io))
        fh.write('#BSUB -P %s\n' %(array_fn))
        fh.write('source %s\n' % self.bashrc)
        fh.write('cd $LS_SUBCWD\n')
        
        # write job array
        for i, fn in enumerate(fns):
            os.chmod(fn, stat.S_IRWXU)
            fh.write('job_array[%d]=%s\n' %(i+1, os.path.abspath(fn)))
        fh.write('${job_array[$LSB_JOBINDEX]};\n')
        fh.close()
        
        # make executable and print message
        os.chmod(array_fn, stat.S_IRWXU)
        message('Writing array %s' %(array_fn))
        return array_fn
    
    def write_PBS_array(self, fns):
        '''write a PBS job array from args and filenames'''
        
        # initialize output file
        fh, array_fn = self.mktemp(suffix='.sh')
        array_fn = os.path.abspath(array_fn)
        
        # write header
        fh.write('#PBS -t 1-%d%%%s\n' %(len(fns), min(len(fns), int(self.l))))
        fh.write('#PBS -e %s.e\n' %(array_fn))
        fh.write('#PBS -q %s\n' %(self.q))
        fh.write('#PBS -o %s.o\n' %(array_fn))
        fh.write('source %s\n' % self.bashrc)
        fh.write('cd $PBS_O_WORKDIR\n')
        
        # write job array
        for i, fn in enumerate(fns):
            os.chmod(fn, stat.S_IRWXU)
            fh.write('job_array[%d]=%s\n' %(i+1, os.path.abspath(fn)))
        fh.write('${job_array[$PBS_ARRAYID]};\n')
        fh.close()
        
        # make executable and print message
        os.chmod(array_fn, stat.S_IRWXU)
        message('Writing array %s' %(array_fn))
        return array_fn
   

    def write_SGE_array(self, fns):
        '''write an SGE job array from args and filenames'''

        # initialize output file
        fh, array_fn = self.mktemp(suffix='.sh')
        array_fn = os.path.abspath(array_fn)

        # write header
        fh.write('#$ -t 1-%d%s\n' %(len(fns), min(len(fns), int(self.l))))
        fh.write('#$ -e %s.e\n' %(array_fn))
        fh.write('#$ -q %s\n' %(self.q))
        fh.write('#$ -o %s.o\n' %(array_fn))
        fh.write('source %s\n' % self.bashrc)
        fh.write('cd $SGE_O_WORKDIR\n')

        # write job array
        for i, fn in enumerate(fns):
            os.chmod(fn, stat.S_IRWXU)
            fh.write('job_array[%d]=%s\n' %(i+1, os.path.abspath(fn)))
        fh.write('${job_array[$SGE_TASK_ID]};\n')
        fh.close()

        # make executable and print message
        os.chmod(array_fn, stat.S_IRWXU)
        message('Writing array %s' %(array_fn))
        return array_fn



 
    def write_job_array(self, fns):
        '''write a job array (LSF or PBS)'''
        if self.cluster == 'broad':
            array_fn = self.write_LSF_array(fns)
	elif self.cluster == 'zcluster':
	    array_fn = self.write_SGE_array(fns)
        elif self.cluster == 'coyote':
            array_fn = self.write_PBS_array(fns)
        return array_fn
    
    def submit(self, commands):
        '''submit a job array to the cluster'''
        fns = self.write_jobs(commands)
        array_fn = self.write_job_array(fns)
        job_ids = self.submit_jobs([array_fn])
        return job_ids
    
    def wait(self, job_ids):
        '''wait for jobs to finish'''
        while True:
            time.sleep(5)
            if self.jobs_finished(job_ids):
                break 
    
    def submit_and_wait(self, commands):
        '''submit job array and wait for it to finish'''
        job_ids = self.submit(commands)
        self.wait(job_ids)
    
    def submit_pipeline(self, pipeline):
        '''a pipeline is a list of lists of commands'''
        for commands in pipeline:
            self.submit_and_wait(commands)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="job submission helper", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('cluster')
    parser.add_argument('queue')
    parser.add_argument('username')
    parser.add_argument('tmp_dir')
    parser.add_argument('bashrc')
    parser.add_argument('--n_cpus', '-n', default=1, type=int, help='number of cpus')
    parser.add_argument('--group', '-G', default=None, help='group')
    parser.add_argument('-m', type=int, default=4, help='memory (gb)')
    parser.add_argument('-l', type=int, default=200, help='job array slot limit')
    parser.add_argument('--io', type=int, default=1, help='disk io (units)')
    parser.add_argument('commands', nargs='?', type=argparse.FileType('r'), default='-')

    # convert the arguments to a dictionary
    args = parser.parse_args()
    argdict = vars(args)
    commands = [line for line in argdict.pop('commands')]

    # create the ssub object with the command line arguments. submit the commands!
    s = Ssub(**argdict)
    s.submit(commands)
