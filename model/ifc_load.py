import numpy as np
import operator

from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepGProp import brepgprop_VolumeProperties, brepgprop_SurfaceProperties, brepgprop_LinearProperties

# http://smartlab1.elis.ugent.be:8889/IFC-repo/
# conda install -c conda-forge -c oce -c dlr-sc -c ifcopenshell ifcopenshell
import ifcopenshell
import ifcopenshell.geom

# RGBA colors for the visualisation of elements
RED, GRAY = (1, 0, 0, 1), (0.6, 0.6, 0.6, 0.1)

# Model freely available at:
# http://www.nibs.org/?page=bsa_commonbimfiles
ifc_file = ifcopenshell.open("../cadfile/Duplex_A_20110907_optimized.ifc")

# Settings to specify usage of pyOCC
settings = ifcopenshell.geom.settings()
settings.set(settings.USE_PYTHON_OPENCASCADE, True)

# Some helper functions to map to the list of walls


def create_shape(elem):
    return ifcopenshell.geom.create_shape(settings, elem)


def calc_volume(s):
    props = GProp_GProps()
    brepgprop_VolumeProperties(s.geometry, props)
    return props.Mass()


def calc_area(s):
    props = GProp_GProps()
    brepgprop_SurfaceProperties(s.geometry, props)
    return props.Mass()


def normalize(li):
    mean = np.mean(li)
    std = np.std(li)
    return map(lambda v: abs(v - mean) / std, li)


# Obtain a list of walls from the model
walls = ifc_file.by_type("IfcWall")
# Create geometry for these walls
shapes = list(map(create_shape, walls))
# Calculate their volumes
volumes = map(calc_volume, shapes)
# Calculate their surface areas
areas = map(calc_area, shapes)
# Compose a feature from the two measures
print(areas)
print(volumes)
feature = normalize(map(operator.floordiv, areas, volumes))

# Initialize the viewer
pyocc_viewer = ifcopenshell.geom.utils.initialize_display()

# Loop over the sorted pairs of feature
# values and corresponding geometry
"""for d, s in sorted(zip(feature, shapes)):
    c = RED if d > 1. else GRAY
    ifcopenshell.geom.utils.display_shape(s, clr=c)
"""
# Fit the model into view
pyocc_viewer.FitAll()

# Allow for user interaction
ifcopenshell.geom.utils.main_loop()
