#!/bin/env python

import sys
import re
import datetime
import subprocess
from collections import defaultdict

RUN_PARAMETERS = "/data/nsc.loki/scripts/lims/getinfo/run-parameters.py"

cols = ['jobid', 'jobname', 'start', 'end', 'comment', 'ncpus']
c = dict((v, i) for i, v in enumerate(cols))

def get_instrument_by_runid(run_id):
    if re.match(r"\d{6}_M", run_id):
        return 'miseq'
    elif re.match(r"\d{6}_NS", run_id):
        return 'nextseq'
    elif re.match(r"\d{6}_[A-Z0-9]", run_id):
        return 'hiseq'
    else:
        return None

sacct_args = [
        "/usr/bin/sacct", "-u", "seq-user", "-P",
        "-S", "0701",  "-o", ",".join(cols), 
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

for row in table:
    if len(row) == len(cols) and row[c['jobname']].endswith("bcl2fastq2"):
        end = datetime.datetime.strptime(row[c['end']], "%Y-%m-%dT%H:%M:%S")
        start = datetime.datetime.strptime(row[c['start']], "%Y-%m-%dT%H:%M:%S")
        time = end - start
        cpus = float(row[c['ncpus']])

        if row[c['comment']] and cpus == 16:
            runinfo = dict(
                    (line.split('|'))
                    for line in 
                    subprocess.check_output(["python", RUN_PARAMETERS, '-p', row[c['comment']]]).splitlines()
                    )
            cycles = [
                    runinfo.get('Read 1 Cycles', "0"),
                    runinfo.get('Read 2 Cycles', "0")
                    ]
            
            instrument = runinfo['Instrument']
            if "Rapid" in runinfo.get("Flow Cell Version", ""):
                instrument += "rapid"

            key = "_".join([instrument] + cycles)
            cpu_hours = time.total_seconds()  / 3600.0
            run_type_cpuhrs[key].append(cpu_hours)
            minutes = time.total_seconds() / 60.0
            run_type_mins[key].append(minutes)

for key in sorted(run_type_cpuhrs.keys()):
    print key, "\t", run_type_cpuhrs[key]

