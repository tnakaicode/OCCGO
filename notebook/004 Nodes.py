#!/usr/bin/env python
# coding: utf-8

# In[22]:


import logging
logging.getLogger("pvtrace").setLevel(logging.CRITICAL)
logging.getLogger("trimesh").setLevel(logging.CRITICAL)
logging.getLogger("matplotlib").setLevel(logging.CRITICAL)
import time
import numpy as np
import functools
import matplotlib
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
from pvtrace import *


# # Nodes
# 
# To position scene objects in useful places you need to `Node` objects.
# 
# They have `translate(dx, dy, dz)` and `rotate(angle, axis)` method which apply to the current pose.
# 
# Let's place the glass box so that the centre is at (1, 1, 1).

# In[23]:


world = Node(
    name="world (air)",
    geometry=Sphere(
        radius=10.0,
        material=Material(refractive_index=1.0)
    )
)
box = Node(
    name="box (glass)",
    geometry=Box(
        size=(1,1,1),
        material=Material(refractive_index=1.5)
    ),
    parent=world
)
box.translate((1, 1, 1))
box.rotate(np.radians(45), (1, 0, 1))
scene = Scene(world)
vis = MeshcatRenderer(wireframe=True)
vis.render(scene)
vis.vis.jupyter_cell()


# Add some widgets for changing the box's location. The sliders change the location tuple of the box node.

# In[24]:


from ipywidgets import interact
import ipywidgets as widgets

def update_box_position(x, y, z):
    box.location = (x, y, z)
    vis.update_transform(box)

interact(
    update_box_position,
    x=widgets.FloatSlider(min=-5,max=5,step=0.1,value=0),
    y=widgets.FloatSlider(min=-5,max=5,step=0.1,value=0),
    z=widgets.FloatSlider(min=-5,max=5,step=0.1,value=0),
)
vis.vis.jupyter_cell()


# ## Set location in scene
# 
# The above sliders are using the box `location` property to set the location of the box in the scene. Here we also need to call the visualisers `update_transform` method to tell it about the change.

# In[25]:


box.location = (-2, 0.0, 0.5)
vis.update_transform(box)  # tell the renderer is need to re-draw


# ## Nested nodes
# 
# An important concept in *pvtrace* is that nodes are nestable and the parent transformation applies defines the childs coordinate system.

# In[27]:


world = Node(
    name="world (air)",
    geometry=Sphere(
        radius=10.0,
        material=Material(refractive_index=1.0)
    )
)
group = Node(
    name="group",
    parent=world
)
box1 = Node(
    name="box 1(glass)",
    geometry=Box(
        size=(1,1,1),
        material=Material(refractive_index=1.5)
    ),
    parent=group
)
box2 = Node(
    name="box 2 (glass)",
    geometry=Box(
        size=(1,1,1),
        material=Material(refractive_index=1.0)
    ),
    parent=group
)
box3 = Node(
    name="box 3 (glass)",
    geometry=Box(
        size=(1,1,1),
        material=Material(refractive_index=1.0)
    ),
    parent=group
)

# Customise location and orientation
box1.location = (-1, 1, 0)
box2.location = (-2, 2, 1)
box3.location = (-3, 3, 2)
group.rotate(np.radians(25), (1, 0, 0))  # applying the rotation to the whole group
scene = Scene(world)
vis = MeshcatRenderer(wireframe=True)
vis.render(scene)
vis.vis.jupyter_cell()


# In[ ]:




