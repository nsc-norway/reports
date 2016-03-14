import os
import sys
import re
import glob
from collections import defaultdict

ignore = [
        ('150821_7001448_0345_BC71MKANXX', 'L001'),
        ]


failed = [
        ('150821_7001448_0345_BC71MKANXX', 'L001'),
        ('151001_D00132_0146_BC7G7GANXX', 'L001'),
        ('151001_D00132_0146_BC7G7GANXX', 'L002'),
        ('151001_D00132_0146_BC7G7GANXX', 'L003'),
        ('151001_D00132_0146_BC7G7GANXX', 'L004')
        ]

inputs = glob.glob("/data/nsc.loki/completed/*/*/Data/Intensities/BaseCalls/QualityControl/Delivery/Email_for_Kent-cowDNA-2015-06-19.xls")
inputs += glob.glob("/data/nsc.loki/issues/*/Data/Intensities/BaseCalls/QualityControl/Delivery/Email_for_Kent-cowDNA-2015-06-19.xls")

sample_reads = defaultdict(int)
sample_reads_pass = defaultdict(int)
total_lanes = set()
total_lanes_pass = set()
sample_lanes_pass = defaultdict(set)


for path in inputs:
    f = open(path)
    run = None
    for line in f:
        columns = line.split("\t")
        if len(columns) == 3 and columns[2].startswith("fragments") and "_R1_" in line:
            if not run:
                print "Error: Run ID not found!"
                sys.exit(1)

            filename = columns[0]
            filename_parts = filename.split("_")
            sample_name = filename_parts[0]
            lane = filename_parts[2]
            
            if (run, lane) in ignore:
                continue

            reads = int(columns[1].replace(",", ""))

            failed_lane = (run, lane) in failed
            print run, "\t", lane, "\t", sample_name, "\t", reads, "\t", "FAILED" if failed_lane else "", "\t"

            sample_reads[sample_name] += reads
            total_lanes.add( (run, lane) )
            if not failed_lane:
                total_lanes_pass.add( (run, lane) )
                sample_lanes_pass[sample_name].add( (run, lane) )
                sample_reads_pass[sample_name] += reads
        else:
            match = re.match(r"Sequence ready for download - sequencing run ([A-Z0-9_]+) - Project_Kent-cowDNA-2015-06-19", line)
            if match:
                run = match.group(1)

print ""
print "-"*30, "SAMPLE TOTALS",  "-"*30
print "" 

print "Sample\tTotal reads\tQC pass reads\tQC pass lanes"
for sample in sorted(sample_reads.keys(), key=lambda sample_name: int(sample_name.split("-")[0])):
    print sample, "\t", sample_reads[sample], "\t", sample_reads_pass[sample], "\t", len(sample_lanes_pass[sample])


print ""
print "-"*30, "GRAND TOTALS",  "-"*30
print "" 
print "Total lanes:       ", len(total_lanes)
print "Total lanes QC OK: ", len(total_lanes_pass)


