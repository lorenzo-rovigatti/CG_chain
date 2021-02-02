#!/usr/bin/env python3

import baggianalysis as ba
import numpy as np
import sys

class CGParser(ba.BaseParser):
    def __init__(self):
        ba.BaseParser.__init__(self)
    
    def _parse_file(self, conf):
        syst = ba.System()
        
        with open(conf) as f:
            # we use the first line to check whether we reached the EOF
            first_line = f.readline().split()
            N = int(first_line[0])
            box = float(first_line[1])
            syst.box = [box, box, box]
            syst.time = int(first_line[2])
            
            if len(first_line) == 0:
                return None
            for _ in range(N):
                x = float(f.readline())
                pos = [x, 0., 0.]
                particle = ba.Particle(syst.available_index(), "0", pos)
                syst.add_particle(particle)
            
        return syst

if len(sys.argv) < 2:
    print("Usage is %s folder" % sys.argv[0], file=sys.stderr)
    exit(1)

parser = CGParser()
trajectory = ba.LazyTrajectory(parser)
trajectory.initialise_from_folder(sys.argv[1], "bead_conf*")

system = trajectory.next_frame()
all_deltas = []
while system != None:
    positions = np.array(system.positions())[:,0]
    delta = positions - np.roll(positions, 1)
    delta[0] += system.box[0]
    all_deltas.append(delta)
    
    system = trajectory.next_frame()

all_deltas = np.concatenate(all_deltas)
hist, bin_edges = np.histogram(all_deltas, bins="auto", density=True)
delta = bin_edges[1] - bin_edges[0]
edges = bin_edges[:-1] + delta / 0.5
histo = np.transpose(np.vstack((edges, hist)))
np.savetxt("PR.dat", histo)
