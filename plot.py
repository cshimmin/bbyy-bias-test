#!/usr/bin/env python

import os
import argparse
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def parse_filename(fname):
    f0 = os.path.basename(fname)
    f0 = f0[:-len('.npy')]
    items = f0.split('-')[1:]
    info = {}
    for itm in items:
        if itm.startswith('x'):
            info['xs'] = float(itm[1:])
        elif itm.startswith('m'):
            info['mass'] = int(itm[1:])
    return info

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--do-fit", action="store_true", help="Do a fit of bias vs. injected")
    parser.add_argument("input_dir", help="The directory containing fit results")
    args = parser.parse_args()

    plt.ion()

    data = {}
    for f in glob(os.path.join(args.input_dir, 'fits-*.npy')):
        info = parse_filename(f)
        xs = info['xs']
        mass = info['mass']

        if xs not in data:
            data[xs] = {}

        if mass not in data[xs]:
            try:
                data[xs][mass] = np.load(f)
            except IOError:
                print "Warning, IOError when opening", f
        else:
            try:
                data[xs][mass] = np.hstack([data[xs][mass], np.load(f)])
            except IOError:
                print "Warning, IOError when opening", f

    adjustments = []
    adjustments_avg = []
    for idx,xs in enumerate(sorted(data.keys(), reverse=True)):
        print "xsec:", xs
        avg = []
        med = []
        masses = sorted(data[xs].keys())
        for mass in masses:
            print "  mass %d  "%mass, len(data[xs][mass])
            avg.append(np.mean(data[xs][mass]))
            med.append(np.median(data[xs][mass]))
        med = np.array(med)
        avg = np.array(avg)

        plt.figure(1)
        plt.plot(masses, med, '.-', color='C%d'%idx, label="inj'd = %g pb"%xs)
        #plt.plot(masses, avg, '--', color='C%d'%idx, label='xs=%g (avg)'%xs)
        plt.figure(2)
        plt.plot(masses, med-xs, '.-', label="inj'd = %g pb (avg: %0.2f +/- %0.2f)"%(xs, np.mean(med-xs), np.std(med-xs)))
        plt.figure(5)
        plt.plot(masses, avg-xs, '.-', label="inj'd = %g pb (avg: %0.2f +/- %0.2f)"%(xs, np.mean(avg-xs), np.std(avg-xs)))
        adjustments.append( (xs, np.mean(med-xs), np.std(med-xs)) )
        adjustments_avg.append( (xs, np.mean(avg-xs), np.std(avg-xs)) )
        if xs>0:
            plt.figure(3)
            plt.plot(masses, (med-xs)/xs, '.-', label="inj'd = %g pb"%xs)
    plt.figure(1)
    plt.ylabel('median best-fit signal [pb]')
    plt.xlabel('mX [GeV]')
    plt.legend()
    plt.savefig(os.path.join(args.input_dir,'median_signal.pdf'))
    plt.figure(2)
    plt.ylabel('median bias (absolute) [pb]')
    plt.xlabel('mX [GeV]')
    plt.legend()
    plt.savefig(os.path.join(args.input_dir,'median_bias.pdf'))
    plt.figure(5)
    plt.ylabel('average bias (absolute) [pb]')
    plt.xlabel('mX [GeV]')
    plt.legend()
    plt.savefig(os.path.join(args.input_dir,'average_bias.pdf'))

    plt.figure(3)
    plt.ylabel('median bias (fractional)')
    plt.xlabel('mX [GeV]')
    plt.legend()
    plt.savefig(os.path.join(args.input_dir,'median_bias_frac.pdf'))
    plt.figure(2)
    plt.figure(1)

    plt.figure()
    adjustments = np.array(sorted(adjustments))
    adjustments_avg = np.array(sorted(adjustments_avg))
    plt.errorbar(adjustments[:,0], adjustments[:,1], yerr=adjustments[:,2], label='medians')
    plt.errorbar(adjustments_avg[:,0], adjustments_avg[:,1], yerr=adjustments_avg[:,2], label='avgs')
    plt.ylabel("avg bias adjust [pb]")
    plt.xlabel("injected xs [pb]")
    plt.legend()

    if args.do_fit:
        def fn(x, a, b, c):
            return a*((1-x)**b) + c
        popt, pcov = curve_fit(fn, adjustments[:,0], adjustments[:,1], sigma=adjustments[:,2])
        print popt
        print pcov
        plt.plot(adjustments[:,0], fn(adjustments[:,0], *popt))

    raw_input('press enter')
