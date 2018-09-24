#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from mpi4py import MPI

comm = MPI.COMM_WORLD
print("Hello! I'm rank %d running on node %s from %d running in total." % (comm.rank, os.uname()[1], comm.size))
comm.Barrier() # wait for sync here