#!/usr/bin/env python
# coding: utf-8

# # Geometry
# 
# This notebook demonstrates the different geometry objects available in pvtrace.
# 
# ## Sphere
# 
# A sphere defined by a radius.
# 
# ## Cylinder
# 
# A cylinder with end caps defined by a radius and a length.
# 
# ## Box
# 
# A box defined by length of the sides.
# 
# ## Mesh
# 
# Any arbitrary mesh can be loaded into pvtrace and traced. However the mesh must be closed for the ray tracing algorithms to work.

# In[3]:


import logging
logging.getLogger('pvtrace').setLevel(logging.CRITICAL)
logging.getLogger('trimesh').setLevel(logging.CRITICAL)
from pvtrace import *
import trimesh


# always need a world node
world = Node(
    name="world",
    geometry=Sphere(
        radius=100.0
    ),
)

# Sphere
sphere = Node(
    name="sphere",
    geometry=Sphere(
        radius=0.5
    ),
    parent=world
)
sphere.translate((0.0, -1.5, 0.5))

# Box
box = Node(
    name="box",
    geometry=Box(
        (1.0, 1.0, 1.0)
    ),
    parent=world
)
box.translate((0.0, 0.0, 0.5))

# Cylinder
cylinder = Node(
    name="cylinder",
    geometry=Cylinder(
        1.0, 0.5
    ),
    parent=world
)
cylinder.translate((0.0, 1.5, 0.5))

# Mesh
mesh = Node(
    name="mesh (icosahedron)",
    geometry=Mesh(
        trimesh=trimesh.creation.icosahedron()
    ),
    parent=world
)
mesh.translate((0.0, 3.5, 1.0))
scene = Scene(world)
vis = MeshcatRenderer(wireframe=True)
vis.render(scene)
vis.vis.jupyter_cell()


# In[ ]:




