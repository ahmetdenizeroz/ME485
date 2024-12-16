import gmsh
import sys
import math

# Before using any functions in the Python API, Gmsh must be initialized:
gmsh.initialize(sys.argv)

# Next add a new model named "cavity" 
gmsh.model.add("mixed")
lc = 0.25

#Points
p1  = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc);

p2  = gmsh.model.geo.addPoint(-0.1,-0.1, 0,lc);
p3  = gmsh.model.geo.addPoint( 0.1,-0.1, 0,lc);
p4  = gmsh.model.geo.addPoint( 0.1, 0.1, 0, lc);
p5  = gmsh.model.geo.addPoint(-0.1, 0.1, 0,lc);

p6  = gmsh.model.geo.addPoint(1/3, 1/3, 0,lc);
p7  = gmsh.model.geo.addPoint(2/3, 2/3, 0,lc);
p8  = gmsh.model.geo.addPoint(-1/3, -1/3, 0,lc);
p9  = gmsh.model.geo.addPoint(-2/3, -2/3, 0,lc);


p10 = gmsh.model.geo.addPoint(-1, -1, 0, lc);
p11 = gmsh.model.geo.addPoint(1, -1, 0, lc);
p12 = gmsh.model.geo.addPoint(1, 1, 0, lc);
p13 = gmsh.model.geo.addPoint(-1, 1, 0, lc);




#As a general rule, elementary entity tags in Gmsh have to be unique per geometrical dimension.
# Inner square
l1  = gmsh.model.geo.addCircleArc(p2, p1, p3)
l2  = gmsh.model.geo.addCircleArc(p3, p1, p4)
l3  = gmsh.model.geo.addCircleArc(p4, p1, p5)
l4  = gmsh.model.geo.addCircleArc(p5, p1, p2)


# Outer boundarties
l5  = gmsh.model.geo.addCircleArc(p10, p1, p11)
l6 = gmsh.model.geo.addCircleArc(p11, p1, p12)
l7 = gmsh.model.geo.addCircleArc(p12, p1, p13)
l8 = gmsh.model.geo.addCircleArc(p13, p1, p10)

l9 = gmsh.model.geo.addCircleArc(p6, p1, p8)
l10 = gmsh.model.geo.addCircleArc(p8, p1, p6)
l11 = gmsh.model.geo.addCircleArc(p7, p1, p9)
l12 = gmsh.model.geo.addCircleArc(p9, p1, p7)




#Surfaces
cl1 = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8]);
cl2 = gmsh.model.geo.addCurveLoop([l11, l12]);
s1 = gmsh.model.geo.addPlaneSurface([cl1,cl2])

c21 = gmsh.model.geo.addCurveLoop([l9, l10]);
c22 = gmsh.model.geo.addCurveLoop([l11, l12]);
s2 = gmsh.model.geo.addPlaneSurface([c21,c22])

c31 = gmsh.model.geo.addCurveLoop([l9, l10]);
c32 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4]);
s3 = gmsh.model.geo.addPlaneSurface([c31,c32])


gmsh.model.geo.mesh.setTransfiniteCurve(l2, 5)
gmsh.model.geo.mesh.setTransfiniteCurve(l3, 5)
gmsh.model.geo.mesh.setTransfiniteCurve(l4, 5)
gmsh.model.geo.mesh.setTransfiniteCurve(l5, 5)


gmsh.model.geo.mesh.setTransfiniteCurve(l9, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l10, 10)

gmsh.model.geo.mesh.setTransfiniteCurve(l5, 25)
gmsh.model.geo.mesh.setTransfiniteCurve(l6, 25)
gmsh.model.geo.mesh.setTransfiniteCurve(l7, 25)
gmsh.model.geo.mesh.setTransfiniteCurve(l8, 25)



"""
# # The `setTransfiniteCurve()' meshing constraints explicitly specifies the
# # location of the nodes on the curve. For example, the following command forces
# # 10 uniformly placed nodes on curve 2 (including the nodes on the two end
# # points):
gmsh.model.geo.mesh.setTransfiniteCurve(l5, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l6, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l7, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l8, 10)

gmsh.model.geo.mesh.setTransfiniteCurve(l1, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l2, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l3, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l4, 10)

gmsh.model.geo.mesh.setTransfiniteCurve(l13, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l14, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l15, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l16, 10)




gmsh.model.geo.mesh.setTransfiniteSurface(s2)
gmsh.model.geo.mesh.setTransfiniteSurface(s3)
gmsh.model.geo.mesh.setTransfiniteSurface(s4)
gmsh.model.geo.mesh.setTransfiniteSurface(s5)




gmsh.model.geo.mesh.setRecombine(2, s2)
gmsh.model.geo.mesh.setRecombine(2, s3)
gmsh.model.geo.mesh.setRecombine(2, s4)
gmsh.model.geo.mesh.setRecombine(2, s5)


"""




gmsh.model.geo.synchronize()


gmsh.model.addPhysicalGroup(1, [l5, l6, l7, l8], 11 , name="outer")
gmsh.model.addPhysicalGroup(1, [l1,l2,l3,l4],       12 , name="inner")
gmsh.model.addPhysicalGroup(2, [s1, s2, s3 ],13 , name="fluid")


# Save it to disk
gmsh.model.mesh.generate(2)
gmsh.write("grad3.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
# This should be called when you are done using the Gmsh Python API:
gmsh.finalize()
