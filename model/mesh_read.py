import meshio

mesh = meshio.read("../cadfile/buckling.ply")
mesh = meshio.read("../cadfile/capacitor.ply")
# mesh = meshio.read("../cadfile/Partition1.ply") # error: need at least one array to concatenate
mesh = meshio.read("../cadfile/tet.ply")
# mesh = meshio.read("../cadfile/unfolding.ply") # error: tuple index out of range
mesh = meshio.read("../cadfile/writeply.ply")
print(mesh)
