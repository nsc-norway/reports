#!/bin/env python

import sys
import re
import datetime
import subprocess
import numpy
from collections import defaultdict

RUN_PARAMETERS = "../run-parameters.py"

cols = ['jobid', 'jobname', 'start', 'end', 'comment', 'ncpus', 'node']
c = dict((v, i) for i, v in enumerate(cols))

sacct_args = [
        "/usr/bin/sacct", "-u", "seq-user", "-P",
        "-S", "0101",  "-o", ",".join(cols), 
        "-s", "COMPLETED"
        ]
try:
    output = subprocess.check_output(sacct_args)
except AttributeError:
    print "Use Python 2.7 for this report"
    print 'Use: scl enable python27 "python demux-report.py"'
    sys.exit(1)


table = [line.split('|') for line in output.split("\n")]

run_type_mins = defaultdict(list)
run_type_cpuhrs = defaultdict(list)

limsid_to_key = {}

run_type_qc_cpuhrs = defaultdict(dict)

for row in table:
    if len(row) == len(cols):
        try:
            end = datetime.datetime.strptime(row[c['end']], "%Y-%m-%dT%H:%M:%S")
            start = datetime.datetime.strptime(row[c['start']], "%Y-%m-%dT%H:%M:%S")
        except ValueError:
            continue
        time = end - start
        cpus = float(row[c['ncpus']])
        
        match = re.match(r"([\d-]+)\.bcl2fastq", row[c['jobname']])
        if match:
            #if cpus != 24.0:
            #    continue

            if row[c['comment']]:# and cpus == 16:
                try:
                    runinfo = dict(
                            (line.split('|'))
                            for line in 
                            subprocess.check_output(["python", RUN_PARAMETERS, '-p', row[c['comment']]]).splitlines()
                            )
                except subprocess.CalledProcessError:
                    continue
                cycles = [
                        runinfo.get('Read 1 Cycles', "0"),
                        runinfo.get('Read 2 Cycles', "0")
                        ]
                
                instrument = runinfo['Instrument']
                if "Rapid" in runinfo.get("Flow Cell Version", ""):
                    instrument += "rapid"

                #key = "_".join([row[c['node']], instrument] + cycles)
                key = "_".join([instrument] + cycles)
                limsid_to_key[match.group(1)] = key
                cpu_hours = cpus * time.total_seconds()  / 3600.0
                run_type_cpuhrs[key].append(cpu_hours)
                minutes = cpus * time.total_seconds() / 60.0
                run_type_mins[key].append(minutes)
        else:
            qc_match = re.match(r"fastqc\.([\d-]+)", row[c['jobname']])
            if qc_match:
                limsid = qc_match.group(1)
                key = limsid_to_key[limsid]
                if not run_type_qc_cpuhrs[key].has_key(limsid):
                    run_type_qc_cpuhrs[key][limsid] = 0
                run_type_qc_cpuhrs[key][limsid] += cpus*time.total_seconds() / 3600.0


print " -- bcl2fastq --"
for key in sorted(run_type_cpuhrs.keys()):
    pad = 20 - len(key)
    print key, " "*pad, "%4.1f +/- %4.1f" % (numpy.mean(run_type_cpuhrs[key]), numpy.std(run_type_cpuhrs[key]))
print " -- fastqc --"
for key in sorted(run_type_qc_cpuhrs.keys()):
    pad = 20 - len(key)
    print key, " "*pad, "%4.1f +/- %4.1f" % (numpy.mean(run_type_qc_cpuhrs[key].values()), numpy.std(run_type_qc_cpuhrs[key].values()))

