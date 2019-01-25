import numpy as np
from unwrap.unwrap import unwrap

from OCC.gp import gp_Pnt, gp_Ax1, gp_Ax3, gp_Vec, gp_Dir, gp_Pnt2d
from OCC.gp import gp_Trsf, gp_Quaternion, gp_Pln, gp_Mat


def pnt_to_xyz(p):
    return p.X(), p.Y(), p.Z()


def occ_to_grasp_rim(axs, pts, filename="pln.rim", name="name", nxy=5):
    pnt = axs.Location()
    trf = gp_Trsf()
    trf.SetTransformation(gp_Ax3(), axs)
    px, py = [], []
    for i in range(len(pts)):
        i0, i1 = i, (i+1) % len(pts)
        p0, p1 = pts[i0].Transformed(trf), pts[i1].Transformed(trf)
        p_x = np.delete(np.linspace(p0.X(), p1.X(), nxy), -1)
        p_y = np.delete(np.linspace(p0.Y(), p1.Y(), nxy), -1)
        px.extend(p_x), py.extend(p_y)

    fp = open(filename, "w")
    fp.write(' {:s}\n'.format(name))
    fp.write('{:12d}{:12d}{:12d}\n'.format(len(px), 1, 1))
    #fp.write(' {:s}\n'.format("mm"))
    for i in range(len(px)):
        data = [px[i], py[i]]
        fp.write(''.join([float_to_string(val) for val in data]) + '\n')


def occ_to_grasp_cor(axs, name="name", filename="pln.cor"):
    pnt = axs.Location()
    v_x = axs.XDirection()
    v_y = axs.YDirection()
    fp = open(filename, "w")
    fp.write(' {:s}\n'.format(name))
    fp.write(' {:s}\n'.format("mm"))
    fp.write(''.join([float_to_string(v) for v in pnt_to_xyz(pnt)]) + '\n')
    fp.write(''.join([float_to_string(v) for v in pnt_to_xyz(v_x)]) + '\n')
    fp.write(''.join([float_to_string(v) for v in pnt_to_xyz(v_y)]) + '\n')
    fp.close()
