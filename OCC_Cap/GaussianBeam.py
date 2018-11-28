import numpy as np
import matplotlib.pyplot as plt

from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Pnt, gp_Vec, gp_Dir, gp_Ax1, gp_Ax2, gp_Ax3
from OCC.gp import gp_Pln, gp_Trsf, gp_Lin, gp_Elips, gp_Elips2d
from OCC.Geom        import Geom_Plane, Geom_Surface, Geom_BSplineSurface
from OCC.Geom        import Geom_Curve, Geom_Line, Geom_Ellipse
from OCC.GeomAPI     import GeomAPI_PointsToBSplineSurface
from OCC.GeomAPI     import GeomAPI_IntCS, GeomAPI_IntSS
from OCC.GeomAPI     import GeomAPI_ProjectPointOnSurf
from OCC.GeomAPI     import GeomAPI_ProjectPointOnCurve
from OCC.GeomAbs     import GeomAbs_C2, GeomAbs_C0, GeomAbs_G1, GeomAbs_G2
from OCC.GeomLProp   import GeomLProp_SurfaceTool
from OCC.GeomProjLib import geomprojlib_Project
from OCC.IntAna      import IntAna_IntConicQuad
from OCC.TopoDS      import TopoDS_HShape, TopoDS_Shape
from OCC.TopoDS      import TopoDS_Face, TopoDS_Wire
from OCC.TopLoc      import TopLoc_Location
from OCC.TColgp      import TColgp_Array2OfPnt
from OCC.BRep           import BRep_Tool
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeFace 
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeWire
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeShell
from OCC.BRepOffsetAPI  import BRepOffsetAPI_ThruSections
from OCC.BRepOffsetAPI  import BRepOffsetAPI_MakePipe
from OCC.BRepPrimAPI    import BRepPrimAPI_MakePrism
from OCCUtils.Construct import make_plane, make_loft
from OCCUtils.Construct import dir_to_vec, vec_to_dir
from OCCUtils.Topology  import Topo, dumpTopology

def surf_trf (axs, face):
    trf = gp_Trsf()
    trf.SetTransformation (axs, gp_Ax3())
    srf = face.Moved (TopLoc_Location (trf))
    return srf

def surf_spl(px, py, pz, axs=gp_Ax3()):
    nx, ny = px.shape
    pnt_2d = TColgp_Array2OfPnt(1, nx, 1, ny)
    for row in range (pnt_2d.LowerRow(), pnt_2d.UpperRow()+1):
        for col in range (pnt_2d.LowerCol(), pnt_2d.UpperCol()+1):
            i,j = row-1, col-1
            pnt = gp_Pnt (px[i, j], py[i, j], pz[i, j])
            pnt_2d.SetValue(row, col, pnt)
    curv = GeomAPI_PointsToBSplineSurface(pnt_2d, 3, 8, GeomAbs_G2, 0.001).Surface()
    surf = BRepBuilderAPI_MakeFace(curv, 1e-6).Face()
    surf = surf_trf (axs, surf)
    return surf

if __name__ == "__main__":
    wave = 1.765
    knum = 2*np.pi / wave
    
    pz = np.linspace (0, 500, 10)
    w0 = 15.0
    rz = pz + 1/pz*(np.pi*w0/wave)**2
    wz = w0 * np.sqrt (1 + (wave*pz/(np.pi*w0**2))**2)
    t0 = wave / (np.pi * w0)
    
    plt.figure()
    plt.plot (pz, rz)

    plt.figure()
    plt.plot (pz, wz)
    plt.plot (pz, np.tan(t0)*pz)
    plt.show ()

    display, start_display, add_menu, add_function_to_menu = init_display()

    api = BRepOffsetAPI_ThruSections()
    for z in np.linspace(0, 1000, 10):
        r_z = z + 1/z*(np.pi*w0/wave)**2
        w_z = w0 * np.sqrt (1 + (wave*z/(np.pi*w0**2))**2)
        print (z, r_z, w_z)
        pnt = gp_Pnt (0, 0, z)
        axs = gp_Ax3 (pnt, gp_Dir(0, 0, 1))
        ax2 = axs.Ax2()
        px  = np.linspace(-1, 1, 100) * 100 
        py  = np.linspace(-1, 1, 100) * 100
        pxy = np.meshgrid (px, py)
        pz  = -1*(pxy[0]**2/(2*r_z) + pxy[1]**2/(2*r_z))
        pln = surf_spl(*pxy, pz, axs)
        wxy = Geom_Ellipse (ax2, w_z, w_z).Elips()
        wxy = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge (wxy).Edge()).Wire() 
        print (wxy, pln)
        api.AddWire(wxy)
        display.DisplayShape (pnt)
        display.DisplayShape (pln)
        display.DisplayShape (wxy)
    api.Build()
    surf_wxy = api.Shape()
    display.DisplayShape(surf_wxy)

    pnt = gp_Pnt (0, 0, 500)    
    axs = gp_Ax3 (pnt, gp_Dir(0, 0, 1))
    axs.Rotate(gp_Ax1(axs.Location(), axs.XDirection()), np.deg2rad(5))
    axs.Rotate(gp_Ax1(axs.Location(), axs.YDirection()), np.deg2rad(30))
    ax2 = axs.Ax2()
    px = np.linspace(-1, 1, 100) * 100 
    py = np.linspace(-1, 1, 100) * 100
    mesh = np.meshgrid (px, py)
    surf = mesh[0]**2/1000 + mesh[1]**2/1000
    #pln = make_plane (axs.Location(), dir_to_vec(axs.Direction()))
    pln = surf_spl(*mesh, surf, axs)
    display.DisplayShape (pln)

    #print (TopoDS_Face(pln))
    print (surf_wxy, pln)
    for face in Topo(surf_wxy).faces():
        surf_wxy = BRep_Tool.Surface(face).GetObject()
    for face in Topo(pln).faces():
        pln = BRep_Tool.Surface(face).GetObject()
    curv = GeomAPI_IntSS (pln.GetHandle(), surf_wxy.GetHandle(), 1.0e-7).Line(1)
    print (curv)
    for p in np.linspace(0, 1, 10):
        curv_p = curv.GetObject().Value(p)
        x = curv_p.X() - gp_Pnt().X()
        y = curv_p.Y() - gp_Pnt().Y()
        z = curv_p.Z() - gp_Pnt().Z()
        r_z = z + 1/z*(np.pi*w0/wave)**2
        w_z = w0 * np.sqrt (1 + (wave*z/(np.pi*w0**2))**2)
        print(x, y, z, r_z, w_z)
        display.DisplayShape (curv_p)
    display.DisplayShape (curv)
    
    display.FitAll()
    start_display()