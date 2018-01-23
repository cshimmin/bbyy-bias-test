#!/usr/bin/env python

import sys, os
import numpy as np
from glob import glob
import cPickle
import numpy as np

if __name__ == "__main__":
    #infile = open(sys.argv[1])
    #times = []
    #for l in infile:
    #    h,m,s = map(int, l.split(":"))
    #    times.append(h*3600 + m*60 + s)
    #print np.mean(times)
    ntotal = 0
    nfail = 0
    nmissing = 0
    input_files = glob(os.path.join(sys.argv[1], "*.pkl"))
    times = []
    processed = []
    for fname in input_files:
        d = cPickle.load(open(fname))
        jid = d['jobid']
        ntotal += len(d['statuses'])
        nfail += max(max(d['statuses']))
        #timeinfo = open('slurm-%s.out'%jid).readlines()[-1]
        #if not timeinfo.startswith("Total time"):
        #    nmissing += 1
        #    continue
        #times.append(float(timeinfo.split(':')[1]))
        times.append(d['runtime'])
        processed.append(len(d['statuses']))
    #print "Failure rate: %.2f%%" % (100.*nfail/ntotal)
    print "Missing time infos:", nmissing
    print 'Total trials:', np.sum(processed)
    print "Avg processed: %.2f (%.2f%%)" % (np.mean(processed), 100.*np.sum(processed)/(len(input_files)*np.max(processed)))
    print "Avg time: %.2f (std=%.2f)"%(np.mean(times), np.std(times))
    print "Median time: %.2f"%np.median(times)
    print "Avg time/trial:", 1.*np.sum(times)/np.sum(processed)
