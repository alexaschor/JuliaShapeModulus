#!/usr/bin/env python3
import sys
import os

# This is a quick and dirty script that will compute the geometry of a DQM Julia set
# in 8x parallel with seq and xargs -P. It relies on trimesh2's mesh_cat (TODO get the
# source for that building in this project). You should adjust the parallel factor
# from 8 to however many cores you have on your machine.

NUM_PARALLEL_JOBS = 8


if len(sys.argv) != 10:
    print("USAGE:")
    print("To create a shaped julia set from a distance field:")
    print(f" {sys.argv[0]} <SDF *.f3d> <P(q) *.poly4d, or the word RANDOM> <output resolution> <a> <b> <offset x> <offset y> <offset z> <output *.obj>")
    print("This will pass all the parameters along to ./bin/run and run it in "+ str(NUM_PARALLEL_JOBS)+"X parallel")
    exit(1)

sdf = sys.argv[1]
p4d = sys.argv[2]
out_res = int(sys.argv[3])
a = sys.argv[4]
b = sys.argv[5]
ox = sys.argv[6]
oy = sys.argv[7]
oz = sys.argv[8]
out_obj = sys.argv[9]

params = [sdf, p4d, str(int(out_res / 2)), a, b, ox, oy, oz, out_obj]

generate_command = "seq 0 7 | xargs -P" + str(NUM_PARALLEL_JOBS) + " -I{} ./bin/run " + " ".join(params) +  ".{}.obj " + "{}"
cat_command = f"mesh_cat {out_obj}.*.obj -o {out_obj}"
rm_command = f"rm {out_obj}.*.obj"

command = "time (" +generate_command + " && " + cat_command + " && " + rm_command + ")"

os.system(command)
