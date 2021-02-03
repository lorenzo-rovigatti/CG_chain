#!/usr/bin/env python3

import baggianalysis as ba
import numpy as np
import sys

class InterBeadFilter(ba.BaseFilter):
    def filter(self, system):
        new_system = system.empty_copy()
        
        positions = np.array(system.positions())[:,0]
        new_positions = positions - np.roll(positions, 1)
        new_positions[0] += system.box[0]
        
        for pos in new_positions:
            new_particle = ba.Particle(new_system.available_index(), "0", [pos, 0., 0.])
            new_system.add_particle(new_particle)
        
        return new_system
        

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
new_filter = InterBeadFilter()
trajectory = ba.LazyTrajectory(parser)
trajectory.add_filter(new_filter)
trajectory.initialise_from_folder(sys.argv[1], "bead_conf*")

msd_obs = ba.MSD(20, False)
msd_obs.analyse_and_print(trajectory, "interbead_msd.dat")
