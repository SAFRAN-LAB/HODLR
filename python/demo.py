#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import time
import kplr
import hodlr
import numpy as np
import matplotlib.pyplot as pl

client = kplr.API()
lc = client.light_curves(12365184, ktc_target_type="SC")[0]
data = lc.read()
t, f, fe, q = (data["TIME"], data["SAP_FLUX"], data["SAP_FLUX_ERR"],
               data["SAP_QUALITY"])
m = np.isfinite(t) * np.isfinite(f) * (q == 0)
t, f, fe, q = t[m], f[m], fe[m], q[m]
mu = np.median(f)
f, fe = f / mu - 1, fe / mu
pl.plot(t, f, ".k", alpha=0.3)
pl.savefig("data.png")

ns = [100, 500, 1000, 5000, 10000, 50000]
builds = np.empty(len(ns))
solves = np.empty(len(ns))
dets = np.empty(len(ns))
for i, n in enumerate(ns):
    print(n)

    # Set up the system.
    diag = fe * fe
    print("Setting up...")
    strt = time.time()
    solver = hodlr.HODLR(1e-8, 0.1, t[:n], diag[:n])
    strt = time.time() - strt
    builds[i] = strt
    print("Took {0} seconds".format(strt))

    # Solve the system.
    print("Solving...")
    strt = time.time()
    xsol = solver.solve(f[:n])
    strt = time.time() - strt
    solves[i] = strt
    print("Took {0} seconds".format(strt))

    # Get the log-determinant.
    print("Finding log-determinant...")
    strt = time.time()
    logdet = solver.logdet()
    strt = time.time() - strt
    dets[i] = strt
    print("Took {0} seconds".format(strt))
    print("Value: {0}".format(logdet))

pl.clf()
pl.plot(np.log10(ns), np.log10(builds+solves+dets), "o")
pl.plot(np.log10(ns), np.log10(np.array(ns) * np.log(ns)), "-")
pl.savefig("builds.png")
