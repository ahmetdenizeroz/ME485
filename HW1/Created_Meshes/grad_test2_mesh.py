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

p6  = gmsh.model.geo.addPoint(-0.5*math.sqrt(2)/2,-0.5*math.sqrt(2)/2, 0,lc);
p7  = gmsh.model.geo.addPoint( 0.5*math.sqrt(2)/2,-0.5*math.sqrt(2)/2, 0,lc);
p8  = gmsh.model.geo.addPoint( 0.5*math.sqrt(2)/2, 0.5*math.sqrt(2)/2, 0,lc);
p9  = gmsh.model.geo.addPoint(-0.5*math.sqrt(2)/2, 0.5*math.sqrt(2)/2, 0,lc);

p10 = gmsh.model.geo.addPoint(0, 0.5, 0,lc);
p11 = gmsh.model.geo.addPoint(-0.5, 0, 0,lc);
p12 = gmsh.model.geo.addPoint(0, -0.5, 0,lc);
p13 = gmsh.model.geo.addPoint(0.5, 0, 0,lc);

p14 = gmsh.model.geo.addPoint(-1, -1, 0, lc);
p15 = gmsh.model.geo.addPoint(1, -1, 0, lc);
p16 = gmsh.model.geo.addPoint(1, 1, 0, lc);
p17 = gmsh.model.geo.addPoint(-1, 1, 0, lc);




#As a general rule, elementary entity tags in Gmsh have to be unique per geometrical dimension.
# Inner square
l1  = gmsh.model.geo.addCircleArc(p2, p1, p3)
l2  = gmsh.model.geo.addCircleArc(p3, p1, p4)
l3  = gmsh.model.geo.addCircleArc(p4, p1, p5)
l4  = gmsh.model.geo.addCircleArc(p5, p1, p2)

# Straight lines to create transfinite part
l5  = gmsh.model.geo.addLine(p6,  p2)
l6  = gmsh.model.geo.addLine(p7 , p3)
l7  = gmsh.model.geo.addLine(p8 , p4)
l8  = gmsh.model.geo.addLine(p9 , p5)
# Outer boundarties
l9  = gmsh.model.geo.addCircleArc(p14, p1, p15)
l10 = gmsh.model.geo.addCircleArc(p15, p1, p16)
l11 = gmsh.model.geo.addCircleArc(p16, p1, p17)
l12 = gmsh.model.geo.addCircleArc(p17, p1, p14)

l13  = gmsh.model.geo.addLine(p8, p10)
l14  = gmsh.model.geo.addLine(p10, p9)
l15  = gmsh.model.geo.addLine(p9, p11)
l16  = gmsh.model.geo.addLine(p11, p6)
l17  = gmsh.model.geo.addLine(p6, p12)
l18  = gmsh.model.geo.addLine(p12, p7)
l19  = gmsh.model.geo.addLine(p7, p13)
l20  = gmsh.model.geo.addLine(p13, p8)


#Surfaces
cl1 = gmsh.model.geo.addCurveLoop([l9, l10, l11, l12]);
cl2 = gmsh.model.geo.addCurveLoop([l13, l14, l15, l16, l17, l18, l19, l20]);
s1 = gmsh.model.geo.addPlaneSurface([cl1,cl2])

cl3 = gmsh.model.geo.addCurveLoop([l13, l14, l8, -l3, -l7]);
s2  = gmsh.model.geo.addPlaneSurface([cl3])

cl4 = gmsh.model.geo.addCurveLoop([-l8, l15, l16, l5, -l4]);
s3  = gmsh.model.geo.addPlaneSurface([cl4])

cl5 = gmsh.model.geo.addCurveLoop([l17, l18, l6, -l1, -l5]);
s4  = gmsh.model.geo.addPlaneSurface([cl5])

cl6 = gmsh.model.geo.addCurveLoop([l19, l20, l7, -l2, -l6]);
s5  = gmsh.model.geo.addPlaneSurface([cl6])


gmsh.model.geo.mesh.setTransfiniteCurve(l3, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l7, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l8, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l14, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l13, 10)

gmsh.model.geo.mesh.setTransfiniteCurve(l1, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l5, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l6, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l18, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l17, 10)


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

"""


gmsh.model.geo.mesh.setRecombine(2, s2)
gmsh.model.geo.mesh.setRecombine(2, s3)
gmsh.model.geo.mesh.setRecombine(2, s4)
gmsh.model.geo.mesh.setRecombine(2, s5)




gmsh.model.geo.synchronize()


gmsh.model.addPhysicalGroup(1, [l9, l10, l11, l12], 11 , name="outer")
gmsh.model.addPhysicalGroup(1, [l1,l2,l3,l4],       12 , name="inner")
gmsh.model.addPhysicalGroup(2, [s1, s2, s3, s4, s5],13 , name="fluid")


# Save it to disk
gmsh.model.mesh.generate(2)
gmsh.write("grad2.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
# This should be called when you are done using the Gmsh Python API:
gmsh.finalize()
