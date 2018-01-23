#!/usr/bin/env python

import sys, os
import argparse
from glob import glob
import cPickle
import matplotlib.pyplot as plt
import numpy as np
A = np.array

def parse_filename(fname):
    f0 = os.path.basename(fname)
    if f0.endswith('.npy.pkl'):
        f0 = f0[:-len('.npy.pkl')]
    else:
        f0 = f0[:-len('.pkl')]
    items = f0.split('-')[1:]
    info = {}
    for itm in items:
        if itm.startswith('x'):
            info['xs'] = float(itm[1:])
        elif itm.startswith('m'):
            info['mass'] = int(itm[1:])
    return info

if __name__ == "__main__":
    print "starting!"
    parser = argparse.ArgumentParser()
    parser.add_argument("--show-all", action="store_true", help="Show all individual pull plots")
    parser.add_argument("--var", default="npbBSM", help="The variable to plot pulls for")
    parser.add_argument("--only-good", action="store_true", help="Only keep trials with 0/0 status")
    parser.add_argument("directory", help='The directory containing pkl files with fit results')
    args = parser.parse_args()

    plt.ion()

    trials = {}
    for fname in glob(os.path.join(args.directory, '*.pkl')):
        info = parse_filename(fname)
        xs = info['xs']
        if xs==0.75: continue
        mass = info['mass']

        if not xs in trials:
            trials[xs] = {}
        if not mass in trials[xs]:
            trials[xs][mass] = []

        pkl = cPickle.load(open(fname))

        status_skip = 0
        for itrial in xrange(len(pkl['vals'])):
            trial_results = {}
            trial_results["statuses"] = pkl['statuses'][itrial]
            if args.only_good and max(trial_results["statuses"])>0:
                status_skip += 1
                continue
            for i,k in enumerate(pkl['keys']):
                trial_results[k] = pkl['vals'][itrial][i]
                trial_results[k+"_lo"] = pkl['errs_lo'][itrial][i]
                trial_results[k+"_hi"] = pkl['errs_hi'][itrial][i]
            trials[xs][mass].append(trial_results)

        if status_skip>0:
            print "Skipped %d trials due to error status. (xs=%g, m=%g)" % (status_skip, xs, mass)

    variable = args.var
    for ixs,xs in enumerate(sorted(trials.keys())):
        masses = A(sorted(trials[xs].keys()))
        pulls = []
        median_bias = []
        for mass in masses:
            median_bias.append(np.median([x['npbBSM']-xs for x in trials[xs][mass]]))
            plt.figure(1)
            plt.clf()
            poi_vals = A([x[variable] for x in trials[xs][mass]])
            if variable=='npbBSM':
                expected=xs
            else:
                expected=0
            poi_errs = np.abs(A([x['%s_hi'%variable] if x[variable]<expected else x['%s_lo'%variable] for x in trials[xs][mass]]))
            poi_pulls = (poi_vals-expected)/poi_errs
            poi_pulls = poi_pulls[poi_errs>0]
            plt.hist(poi_pulls, histtype='step', label='avg=%0.2f med=%0.2f std=%0.2f'%(np.mean(poi_pulls), np.median(poi_pulls), np.std(poi_pulls)))
            pulls.append( (np.mean(poi_pulls), np.std(poi_pulls)) )
            plt.axvline(0, color='black')
            plt.title("%s (xs=%g, m=%g)" % (variable, xs, mass))
            plt.legend()
            plt.figure(2)
            plt.clf()
            plt.hist(poi_vals, histtype='step', label='avg=%0.2f med=%0.2f std=%0.2f'%(np.mean(poi_vals), np.median(poi_pulls), np.std(poi_vals)))
            plt.hist(poi_vals[poi_errs>0], histtype='step', label='(nonzero errro) avg=%0.2f std=%0.2f'%(np.mean(poi_vals[poi_errs>0]), np.std(poi_vals[poi_errs>0])))
            plt.title("%s vals (xs=%g, m=%g)" % (variable, xs, mass))
            plt.legend()
            plt.figure(3)
            plt.clf()
            plt.hist(poi_errs, histtype='step', label='avg=%0.2f std=%0.2f'%(np.mean(poi_errs), np.std(poi_errs)))
            plt.title("%s errs (xs=%g, m=%g)" % (variable, xs, mass))
            plt.legend()
            if args.show_all:
                raw_input("press enter")
        pulls = A(pulls)

        plt.figure(4)
        if ixs==0:
            plt.fill_between(masses, 1.*np.ones_like(masses), -1.*np.ones_like(masses), facecolor='black', alpha=0.2)
            plt.plot(masses, np.zeros_like(masses), color='black')
        plt.errorbar(masses+2*ixs, pulls[:,0], yerr=pulls[:,1], fmt='o', label='xs=%g'%xs)
        plt.title("%s pulls"%variable)
        plt.legend(loc='lower center')
        
        plt.figure(5)
        plt.plot(masses, median_bias, label='xs=%g'%xs)
        plt.title('median bias')
        plt.legend()
        raw_input('press enter')
    raw_input('press enter')
