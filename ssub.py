#!/usr/bin/env python

import argparse, os, re, select, stat, subprocess, sys, tempfile, time
from util import *

'''
ssub is a simple job submission script

usage:
  cat list_of_commands.txt | ssub -n 100 -q hour -G gscidfolk -m 8 --io 10
  ssub -n 100 -q week -G broadfolk -m 8 --io 10 "command1; command2; command3;"

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

# set global variables
config = ConfigParser.ConfigParser()
config.read('train.cfg')
username = config.get('User', 'username')
temp_dir = config.get('User', 'tmp_directory')
cluster= config.get('User', 'cluster')


def parse_args():
    # print usage statement
    usage = "cat list_of_commands.txt | ssub -n 100 -q hour -G gscidfolk -m 8 --io 10\n"
    usage +="ssub -n 100 -q week -G broadfolk -m 8 --io 10 'command 1; command 2; command 3;"
    
    # add command line arguments
    parser = argparse.ArgumentParser(usage = usage)
    parser.add_argument('-n', default = -1, type = int, help = 'number of cpus')
    parser.add_argument('-q', default = 'hour', help = 'queue')
    parser.add_argument('-G', default = 'gscidfolk', help = 'group')
    parser.add_argument('-m', default = 4, help = 'memory (gb)')
    parser.add_argument('-l', default = 200, help = 'job array slot limit')
    parser.add_argument('--io', default = 1, help = 'disk io (units)')
    parser.add_argument('commands', nargs = '?', default = '')
    
    # parse arguments from stdin
    if __name__ == '__main__':
        args = parser.parse_args()
    else:
        args = parser.parse_args('')
    return args


class Ssub():
    
    def __init__(self, cluster = 'coyote'):
        
        # get command line arguments
        args = parse_args()
        
        # initialize cluster parameters
        self.cluster = cluster
        self.username = username
        self.temp_dir = temp_dir
        self.header = '#!/bin/bash\n'
        self.l = args.l
        self.n = args.n
        self.m = args.m
        self.q = args.q
        self.G = args.G
        self.io = args.io
        self.commands = args.commands
        
        # broad parameters
        if cluster == 'broad':
            self.submit_cmd = 'bsub'
            self.stat_cmd = 'bjobs -w'
            self.parse_job = lambda x: re.search('Job <(\d+)>', x).group(1)
        
        # coyote parameters
        elif cluster == 'coyote':
            self.submit_cmd = 'qsub'
            self.stat_cmd = 'qstat -l'
            self.parse_job = lambda x: x.rstrip()
        
        # unrecognized cluster
        else:
            error('unrecognized cluster %s' %(cluster))
    
            
    def __repr__(self):
        return '\ncluster: %s\nusername: %s\ntemp_dir: %s\nn: %d\nm: %d\nq: %s\nG: %s\nio: %s\n' %(self.cluster, self.username, self.temp_dir, self.n, self.m, self.q, self.G, self.io)
    
    
    def __str__(self):
        a = self.__repr__()
        a = re.sub('\n', '\n  ', a)
        return a
    
    
    def mktemp(self, suffix = '.tmp'):
        # make temporary file and return [fh, fn]
        fh, fn = tempfile.mkstemp(suffix='.sh', dir=self.temp_dir)
        os.close(fh)
        fh = open(fn, 'w')
        fh.write(self.header)
        return fh, fn
    
    
    def job_status(self):
        # return job status
        process = subprocess.Popen([self.stat_cmd], stdout = subprocess.PIPE, shell=True)
        [out, err] = process.communicate()
        out = [line for line in out.split('\n') if self.username in line]
        return [out, err]
    
    
    def jobs_finished(self, job_ids):
        # check if jobs are finished
        [out, err] = self.job_status()
        for job in out:
            for job_id in job_ids:
                if job_id in job:
                    return False
        return True
    
    
    def write_jobs(self, commands):
        # write job scripts from a list of commands
        
        # initialize output files
        fhs, fns = zip(*[self.mktemp(suffix='.sh') for i in range(min(self.n, len(commands)))])
        fhs_cycle = cycle(fhs)
        
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
            message('Writing job %s' %(fn))
        
        return fns
    
    
    def submit_jobs(self, fns):
        # submit jobs to the cluster
        job_ids = []
        for fn in fns:
            if self.cluster == 'broad':
                process = subprocess.Popen(['%s < %s' %(self.submit_cmd, fn)], stdout = subprocess.PIPE, shell=True)
            elif self.cluster == 'coyote':
                process = subprocess.Popen(['%s %s' %(self.submit_cmd, fn)], stdout = subprocess.PIPE, shell=True)
            else:
                quit()
            [out, error] = process.communicate()
            job_ids.append(self.parse_job(out))
            message('Submitting job %s' %(fn))
        return job_ids
    
    
    def write_LSF_array(self, fns):
        # write an LSF job array from args and filenames
        
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
        fh.write('source /home/unix/%s/.bashrc\n' %(self.username))
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
        # write a PBS job array from args and filenames
        
        # initialize output file
        fh, array_fn = self.mktemp(suffix='.sh')
        array_fn = os.path.abspath(array_fn)
        
        # write header
        fh.write('#PBS -t 1-%d%%%s\n' %(len(fns), min(len(fns), int(self.l))))
        fh.write('#PBS -e %s.e\n' %(array_fn))
        fh.write('#PBS -o %s.o\n' %(array_fn))
        fh.write('source /home/%s/.bashrc\n' %(self.username))
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
    
    
    def write_job_array(self, fns):
        # write a job array (LSF or PBS)
        if self.cluster == 'broad':
            array_fn = self.write_LSF_array(fns)
        elif self.cluster == 'coyote':
            array_fn = self.write_PBS_array(fns)
        return array_fn
    
    
    def submit(self, commands, out = False):
        # submit a job array to the cluster
        if out == False:
            fns = self.write_jobs(commands)
            array_fn = self.write_job_array(fns)
            job_ids = self.submit_jobs([array_fn])
            return job_ids
        elif out == True:
            print '\n'.join(commands)
            return []
    
    
    def wait(self, job_ids, out = False):
        # wait for jobs to finish
        if out == False:
            while True:
                time.sleep(5)
                if self.jobs_finished(job_ids):
                    break 
    
    
    def submit_and_wait(self, commands, out = False):
        # submit job array and wait for it to finish
        job_ids = self.submit(commands, out = out)
        self.wait(job_ids, out = out)
    
    
    def submit_pipeline(self, pipeline, out = False):
        # a pipeline is a list of lists of commands
        for commands in pipeline:
            self.submit_and_wait(commands, out = out)
    
    
    def validate_output(self, fns, out = False):
        # make sure files exist
        if out == False:
            test = [os.path.exists(fn) for fn in fns]
            if False in test:
                error('file %s does not exist' %(fns[test.index(False)]))
    
    
    def remove_files(self, fns, out = False, to_queue = False):
        # remove list of files
        cmds = ['rm %s' %(fn) for fn in fns]
        if out == False:
            if to_queue == False:
                for cmd in cmds:
                    os.system(cmd)
            elif to_queue == True:
                self.submit_and_wait(cmds, out = out)
        elif out == True:
            print '\n'.join(cmds)
    
    
    def move_files(self, x, y, out = False, to_queue = False):
        # move files from x to y
        if len(x) != len(y):
            error('move_files: len(x) != len(y)')
        cmds = ['mv %s %s' %(a, b) for a, b in zip(x, y)]
        if out == False:
            if to_queue == False:
                for cmd in cmds:
                    os.system(cmd)
            elif to_queue == True:
                self.submit_and_wait(cmds, out = out)
        else:
            print '\n'.join(cmds)
    
    
    def gzip_and_validate(self, fns, tgz, out = False, to_queue = False):
        # compress files and validate .gz file
        cmd = 'tar -cvzf %s %s' %(tgz, ' '.join(fns))
        if out == False:
            if to_queue == False:
                os.system(cmd)
            else:
                self.submit_and_wait(args, [cmd])
            try:
                test = subprocess.check_output(['gunzip', '-t', tgz])
                if test != '':
                    error('gunzip -t %s failed' %(tgz))
            except:
                error('gunzip -t %s failed' %(tgz))
        elif out == True:
            print cmd
    
    
    def run_local(self, commands, out = False):
        # Run commands locally
        if out == False:
            for command in commands:
                os.system(command)
        elif out == True:
            print '\n'.join(commands)
    


def initialize():
    # initialize global variables for ssub
    
    # parse command line args
    ssub = Ssub()
    
    # get list of commands
    commands = []
    if ssub.commands != '':
        commands += [command.strip() for command in ssub.commands.split(';')]
    if select.select([sys.stdin], [], [], 0)[0]:
        commands += [line.rstrip() for line in sys.stdin.readlines()]
    
    # calculate number of cpus
    if ssub.n < 0:
        ssub.n = len(commands)
    
    return ssub, commands


ssub, commands = initialize()

if __name__ == '__main__':
    ssub.submit(commands)
