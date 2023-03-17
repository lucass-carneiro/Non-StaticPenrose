#!/usr/bin/python

import os
import subprocess

def make_file_list(start, end):
    list = []
    for i in range(start, end + 1):
        list.append("sks_breakup_M2_" +  str(i).zfill(3))
    return list

def run_batch(start, end):
    files = make_file_list(start, end)
    proc_data = []
    
    for file in files:
        os.mkdir(file)
        os.chdir(file)

        os.symlink("../Release/bin/grlensing", "grlensing")
        os.symlink("../Release/lib/", "lib")
        os.symlink("../configs/grlensing_config.yaml", "grlensing_config.yaml")
        os.symlink("../configs/sks/" + file + ".yaml", file + ".yaml")

        stdout_file = open("stdout.txt", "w")
        stderr_file = open("stderr.txt", "w")
        
        proc = subprocess.Popen([
            "mpirun",
            "-n",
            "3",
            "grlensing",
            "penrose-breakup",
            file + ".yaml"
            ],
            stdout=stdout_file,
            stderr=stderr_file
        )
        
        proc_data.append([proc, stdout_file, stderr_file])

        os.chdir("../")
    
    for pd in proc_data:
         pd[0].communicate()
         pd[1].close()
         pd[2].close()

run_batch(18, 23)
run_batch(24, 29)
run_batch(30, 35)
run_batch(36, 41)
run_batch(42, 47)
run_batch(48, 53)
run_batch(54, 59)