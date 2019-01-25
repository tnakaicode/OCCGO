echo off

python surf_coord_abs.py --surf=surf1_beam --pxyz=0 0 0 --rxyz=0.1 0.1 0

python surf_coord_abs.py --surf=surf1 --pxyz=0 0 0 --rxyz=0 0 0
python surf_coord_abs.py --surf=surf2 --pxyz=0 0 100 --rxyz=0 0 0
python surf_coord_abs.py --surf=surf3 --pxyz=0 100 100 --rxyz=0 0 0
python surf_coord_abs.py --surf=surf4 --pxyz=0 100 200 --rxyz=0 0 0

python surf.py --surf=surf1 --lxy=100 100 --nxy=200 200 --rxy=0 0
python surf.py --surf=surf2 --lxy=100 100 --nxy=200 200 --rxy=100 0
python surf.py --surf=surf3 --lxy=100 100 --nxy=200 200 --rxy=100 0
python surf.py --surf=surf4 --lxy=100 100 --nxy=200 200 --rxy=0 0

python surf_stp.py --surf=surf1
python surf_stp.py --surf=surf2
python surf_stp.py --surf=surf3
python surf_stp.py --surf=surf4
