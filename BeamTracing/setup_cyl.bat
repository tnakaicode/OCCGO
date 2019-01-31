echo off

python surf_coord_abs.py --surf=beam_cyl --pxyz=9.0 0.0 0.0 --rxyz=90.0 0.0 0.0
python surf_coord_rel.py --surf=beam_cyl --pxyz=0.0 0.0 0.0 --rxyz=0.0 -90.0 0.0
python surf_coord_rel.py --surf=beam_cyl --pxyz=0.0 0.0 0.0 --rxyz=-0.1 75.0 0.0
python surf_cylinder.py --lxy=-1.0 180.0 --radi=10.0 10.0 --rxy=2.0 2.0
python ray_trace_cyl.py