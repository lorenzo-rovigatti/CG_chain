#!/usr/bin/env python3

import baggianalysis as ba

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

parser = CGParser()
trajectory = ba.LazyTrajectory(parser)
trajectory.initialise_from_folder("confs", "bead_conf*")

msd_obs = ba.MSD(20, True)
msd_obs.analyse_and_print(trajectory, "msd.dat")
