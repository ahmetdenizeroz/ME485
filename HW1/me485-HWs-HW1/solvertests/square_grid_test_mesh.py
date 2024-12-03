import gmsh
import sys

gmsh.initialize(sys.argv)
gmsh.model.add("squaregrid")

# Parameters
region_size = 3.0  # Size of the entire square region
num_divisions = 3  # Number of divisions along each edge

# Create points for the large square region
p1 = gmsh.model.geo.addPoint(0, 0, 0)
p2 = gmsh.model.geo.addPoint(region_size, 0, 0)
p3 = gmsh.model.geo.addPoint(region_size, region_size, 0)
p4 = gmsh.model.geo.addPoint(0, region_size, 0)

# Create lines for the large square region
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)

# Create curve loop and surface for the large square region
cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
surface = gmsh.model.geo.addPlaneSurface([cl])

# Apply transfinite meshing to lines and surface
gmsh.model.geo.mesh.setTransfiniteCurve(l1, num_divisions + 1)
gmsh.model.geo.mesh.setTransfiniteCurve(l2, num_divisions + 1)
gmsh.model.geo.mesh.setTransfiniteCurve(l3, num_divisions + 1)
gmsh.model.geo.mesh.setTransfiniteCurve(l4, num_divisions + 1)
gmsh.model.geo.mesh.setTransfiniteSurface(surface, "Left")

gmsh.model.geo.mesh.setRecombine(2, surface)

# Synchronize the model
gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(1, [l1, l2, l3, l4], 11, name = "outer")
gmsh.model.addPhysicalGroup(2, [surface], 13, name = "fluid")

# Set the meshing algorithm to structured
#gmsh.option.setNumber("Mesh.Algorithm", 5)  # Use the transfinite algorithm

# Generate 2D mesh
gmsh.model.mesh.generate(2)

# Save the mesh
gmsh.write("squaregrid.msh")

# Optionally display the GUI
gmsh.fltk.run()

# Finalize Gmsh
gmsh.finalize()
