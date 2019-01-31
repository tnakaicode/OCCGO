echo off

python surf_coord_abs.py --surf=beam_cyl --pxyz=9.0 0.0 5.0 --rxyz=90.0 0.0 0.0
python surf_coord_rel.py --surf=beam_cyl --pxyz=0.0 0.0 0.0 --rxyz=0.0 -90.0 0.0
python surf_coord_rel.py --surf=beam_cyl --pxyz=0.0 0.0 0.0 --rxyz=-1.0 100.0 0.0
python surf_cylinder.py --lxy=0.0 200.0
python ray_trace_cyl.py