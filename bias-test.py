#!/usr/bin/env python

import time
START_TIME = time.time()

import argparse
import ROOT as r
r.gROOT.SetBatch(1)
import sys, os
import numpy as np
import cPickle


def iterset(rooset):
    itr = rooset.createIterator()
    while True:
        x = itr.Next()
        if not x: break
        yield x

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ws", required=True, help="The workspace file")
    parser.add_argument("--seed", type=int, default=1, help="The random seed for RooFit")
    parser.add_argument("--out", help="Output filename")
    parser.add_argument("--ntrial", type=int, default=10, help="The number of trials generate")
    parser.add_argument("--mX", type=int, default=300, help="The resonance mass [GeV]")
    parser.add_argument("--xsec", type=float, default=1.0, help="The signal cross section to inject [pb]")
    parser.add_argument("--freeze-bias", action="store_true", help="Fix the BIAS NP at zero")
    parser.add_argument("--free-shape", action="store_true", help="Use free-floating shape params")
    parser.add_argument("--freeze-shape", action="store_true", help="Fix the shape params to constant values")
    parser.add_argument("--free-norm", action="store_true", help="Use free-floating norm params")
    parser.add_argument("--bias-adj", type=float, help="Apply a bias adjust offset")
    parser.add_argument("--poisson", action="store_true", help="Randomize number of generated events by poisson sampling.")
    parser.add_argument("--only-good", action="store_true", help="Only write out fits that had 0/0 status.")
    parser.add_argument("--skip-minos", action="store_true", help="Do not run minos, only migrad")
    parser.add_argument("--hesse", action="store_true", help="Run Hesse after Migrad")
    args = parser.parse_args()

    r.RooRandom.randomGenerator().SetSeed(args.seed)
    np.random.seed(args.seed+100)

    # WARNING ERROR FATAL
    r.RooMsgService.instance().setGlobalKillBelow(r.RooFit.ERROR)

    f = r.TFile(args.ws)
    w = f.Get("combination")
    mc = w.obj("mconfig")

    pdf = mc.GetPdf()

    if not args.bias_adj is None:
        w.factory("expr::npbBSM_adj('npbBSM-bias_adj',npbBSM,bias_adj[0])")
        w.factory("EDIT::pdf_alt(%s,npbBSM=npbBSM_adj)"%pdf.GetName())
        w.obj("bias_adj").setVal(args.bias_adj)
        pdf = w.obj("pdf_alt")

    obs = w.obj("gg_mass")
    cat = w.obj("channellist")

    mX = w.obj("mHiggs")
    xsec = w.obj("npbBSM")

    # set NPs to zero
    for v in iterset(mc.GetNuisanceParameters()):
        print "setting %s=0"%v.GetName()
        v.setVal(0)

    if args.freeze_bias:
        print "fixing BIAS=0"
        w.obj("BIAS").setVal(0)
        w.obj("BIAS").setConstant(True)

    if args.free_shape or args.freeze_shape:
        if args.free_shape:
            print "setting shape parameters to free"
            w.obj("novosibirsk_peak_bj").setConstant(False)
            w.obj("novosibirsk_peak_bb").setConstant(False)
            w.obj("novosibirsk_tail_bj").setConstant(False)
            w.obj("novosibirsk_tail_bb").setConstant(False)
            w.obj("novosibirsk_width_bj").setConstant(False)
            w.obj("novosibirsk_width_bb").setConstant(False)
        else:
            print "Fixing shape parameters to constant"
        w.obj("bkg_constraint_shape_bj").setConstant(True)
        w.obj("bkg_constraint_shape_bj").setVal(0)
        w.obj("bkg_constraint_shape_bb").setConstant(True)
        w.obj("bkg_constraint_shape_bb").setVal(0)
        w.obj("bkg_constraint_tail_bj").setConstant(True)
        w.obj("bkg_constraint_tail_bj").setVal(0)
        w.obj("bkg_constraint_tail_bb").setConstant(True)
        w.obj("bkg_constraint_tail_bb").setVal(0)
        w.obj("bkg_constraint_width_bj").setConstant(True)
        w.obj("bkg_constraint_width_bj").setVal(0)
        w.obj("bkg_constraint_width_bb").setConstant(True)
        w.obj("bkg_constraint_width_bb").setVal(0)

    if args.free_norm:
        w.obj("bkg_constraint_bj").setConstant(True)
        w.obj("bkg_constraint_bj").setVal(0)
        w.obj("bkg_constraint_bb").setConstant(True)
        w.obj("bkg_constraint_bb").setVal(0)
        w.obj("nbkg_fit_bj_bj").setConstant(False)
        w.obj("nbkg_fit_bb_bb").setConstant(False)

    print "setting mX=%g"%args.mX
    mX.setVal(args.mX)
    print "fixing xsec=%g"%args.xsec
    xsec.setVal(args.xsec)
    xsec.setConstant(True)

    # generate datasets
    datasets = []
    expected_events = pdf.expectedEvents(r.RooArgSet(cat,obs))
    for itrial in xrange(args.ntrial):
        if args.poisson:
            ds = pdf.generate(r.RooArgSet(cat, obs), np.random.poisson(expected_events))
        else:
            ds = pdf.generate(r.RooArgSet(cat, obs))
        ds.SetName("ds_%03d"%itrial)
        datasets.append(ds)

    xsec.setVal(1)
    xsec.setConstant(False)

    other_keys = ['npbBSM', 'bias_bj', 'bias_bb']
    if args.free_norm:
        other_keys.extend([
            'nbkg_fit_bj_bj',
            'nbkg_fit_bb_bb',
            ])
    else:
        other_keys.extend([
            'bkg_constraint_bj',
            'bkg_constraint_bb',
            ])
    if args.free_shape:
        other_keys.extend([
            'novosibirsk_peak_bj',
            'novosibirsk_peak_bb',
            'novosibirsk_tail_bj',
            'novosibirsk_tail_bb',
            'novosibirsk_width_bj',
            'novosibirsk_width_bb',
            ])
    elif not args.freeze_shape:
        other_keys.extend([
            'bkg_constraint_shape_bj',
            'bkg_constraint_shape_bb',
            'bkg_constraint_tail_bj',
            'bkg_constraint_tail_bb',
            'bkg_constraint_width_bj',
            'bkg_constraint_width_bb',
            ])

    poi_vals = []
    other_vals = []
    errs_lo = []
    errs_hi = []
    statuses = []
    status_skip = 0
    for ds in datasets:
        nll = pdf.createNLL(ds)

        fit_statuses = []
        minimizer = r.RooMinuit(nll)
        minimizer.migrad()
        fit_statuses.append(minimizer.save().status())
        if args.hesse:
            minimizer.hesse()
            fit_statuses.append(minimizer.save().status())
        if not args.skip_minos:
            minimizer.minos()
            fit_statuses.append(minimizer.save().status())

        if args.only_good and max(fit_statuses)>0:
            status_skip += 1
            continue

        statuses.append( tuple(fit_statuses) )

        poi_vals.append(xsec.getVal())
        other_vals.append(tuple([w.obj(k).getVal() for k in other_keys]))
        errs_lo.append(tuple([w.obj(k).getAsymErrorLo() for k in other_keys]))
        errs_hi.append(tuple([w.obj(k).getAsymErrorHi() for k in other_keys]))
        if args.out:
            np.save(args.out, poi_vals)
            with open(args.out+'.pkl', 'wb') as fpkl:
                cPickle.dump(dict(
                        keys=other_keys,
                        vals=other_vals,
                        errs_lo=errs_lo,
                        errs_hi=errs_hi,
                        statuses=statuses,
                        argv=sys.argv,
                        jobid=os.environ.get('SLURM_JOBID', None),
                        runtime=(time.time() - START_TIME)
                    ),
                    fpkl)

    print "Skipped %d trials with bad status."%status_skip
    print "Total time:", (time.time() - START_TIME)
