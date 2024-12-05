import gmsh
import math
import sys

# Initialize Gmsh
gmsh.initialize()
gmsh.model.add("mixed")
lc = 0.25  # Characteristic length for meshing

# Points
# Center point
p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc)

# Inner octagonal
p2 = gmsh.model.geo.addPoint(0.4, 0.4 + 0.8 * math.cos(math.pi/6), 0.0, lc)
p3 = gmsh.model.geo.addPoint(-0.4, 0.4 + 0.8 * math.cos(math.pi/6), 0.0, lc)
p4 = gmsh.model.geo.addPoint(-0.8, 0.4 , 0.0, lc)
p5 = gmsh.model.geo.addPoint(-0.8, -0.4 , 0.0, lc)
p6 = gmsh.model.geo.addPoint(-0.4, -0.4 - 0.8 * math.cos(math.pi/6), 0.0, lc)
p7 = gmsh.model.geo.addPoint(0.4, -0.4 - 0.8 * math.cos(math.pi/6), 0.0, lc)
p8 = gmsh.model.geo.addPoint(0.8, -0.4 , 0.0, lc)
p9 = gmsh.model.geo.addPoint(0.8, 0.4 , 0.0, lc)

# Inner circle points   
p10 = gmsh.model.geo.addPoint(0.4/2, (0.4 + 0.8 * math.cos(math.pi/6))/2, 0.0, lc)
p11 = gmsh.model.geo.addPoint(-0.4/2, (0.4 + 0.8 * math.cos(math.pi/6))/2, 0.0, lc)
p12 = gmsh.model.geo.addPoint(-0.8/2, 0.4/2, 0.0, lc)
p13 = gmsh.model.geo.addPoint(-0.8/2, -0.4/2, 0.0, lc)
p14 = gmsh.model.geo.addPoint(-0.4/2, (-0.4 - 0.8 * math.cos(math.pi/6))/2, 0.0, lc)
p15 = gmsh.model.geo.addPoint(0.4/2, (-0.4 - 0.8 * math.cos(math.pi/6))/2, 0.0, lc)
p16 = gmsh.model.geo.addPoint(0.8/2, -0.4/2, 0.0, lc)
p17 = gmsh.model.geo.addPoint(0.8/2, 0.4/2, 0.0, lc)

# Outer circle points (forming a proper circle)
p18 = gmsh.model.geo.addPoint(2.0, 2.0, 0.0, lc)
p19 = gmsh.model.geo.addPoint(-2.0, 2.0, 0.0, lc)
p20 = gmsh.model.geo.addPoint(-2.0, -2.0, 0.0, lc)
p21 = gmsh.model.geo.addPoint(2.0, -2.0, 0.0, lc)

# Lines connecting octagonal vertices to inner circle points
l1 = gmsh.model.geo.addLine(p2, p10)
l2 = gmsh.model.geo.addLine(p3, p11)
l3 = gmsh.model.geo.addLine(p4, p12)
l4 = gmsh.model.geo.addLine(p5, p13)
l5 = gmsh.model.geo.addLine(p6, p14)
l6 = gmsh.model.geo.addLine(p7, p15)
l7 = gmsh.model.geo.addLine(p8, p16)
l8 = gmsh.model.geo.addLine(p9, p17)


# Lines connecting octagonal vertices 
l9 = gmsh.model.geo.addLine(p2, p3)
l10 = gmsh.model.geo.addLine(p3, p4)
l11 = gmsh.model.geo.addLine(p4, p5)
l12 = gmsh.model.geo.addLine(p5, p6)
l13 = gmsh.model.geo.addLine(p6, p7)
l14 = gmsh.model.geo.addLine(p7, p8)
l15 = gmsh.model.geo.addLine(p8, p9)
l16 = gmsh.model.geo.addLine(p9, p2)


# Inner hexagonal points ( forgive my typo, miscalculation)
arc1 = gmsh.model.geo.addLine(p10, p11)
arc2 = gmsh.model.geo.addLine(p11, p12)
arc3 = gmsh.model.geo.addLine(p12, p13)
arc4 = gmsh.model.geo.addLine(p13, p14)
arc5 = gmsh.model.geo.addLine(p14, p15)
arc6 = gmsh.model.geo.addLine(p15, p16)
arc7 = gmsh.model.geo.addLine(p16, p17)
arc8 = gmsh.model.geo.addLine(p17, p10)

# Outer circle arcs
arc9 = gmsh.model.geo.addCircleArc(p18, p1, p19)
arc10 = gmsh.model.geo.addCircleArc(p19, p1, p20)
arc11 = gmsh.model.geo.addCircleArc(p20, p1, p21)
arc12 = gmsh.model.geo.addCircleArc(p21, p1, p18)


# Create curve loops for surfaces
cl1 = gmsh.model.geo.addCurveLoop([l9,l10,l11,l12,l13,l14,l15,l16])
cl2 = gmsh.model.geo.addCurveLoop([arc9, arc10, arc11, arc12])
s1 = gmsh.model.geo.addPlaneSurface([cl1,cl2])


cl3 = gmsh.model.geo.addCurveLoop([l9,l2,-arc1,-l1])
cl4= gmsh.model.geo.addCurveLoop([l10,l3,-arc2,-l2])
cl5 = gmsh.model.geo.addCurveLoop([l11, l4, -arc3, -l3])
cl6 = gmsh.model.geo.addCurveLoop([l12, l5, -arc4, -l4])
cl7 = gmsh.model.geo.addCurveLoop([l13, l6, -arc5, -l5])
cl8 = gmsh.model.geo.addCurveLoop([l14, l7, -arc6, -l6])
cl9 = gmsh.model.geo.addCurveLoop([l15, l8, -arc7, -l7])
cl10 = gmsh.model.geo.addCurveLoop([l16, l1, -arc8, -l8])

s2 = gmsh.model.geo.addPlaneSurface([cl3])
s3 = gmsh.model.geo.addPlaneSurface([cl4])
s4 = gmsh.model.geo.addPlaneSurface([cl5])
s5 = gmsh.model.geo.addPlaneSurface([cl6])
s6 = gmsh.model.geo.addPlaneSurface([cl7])
s7 = gmsh.model.geo.addPlaneSurface([cl8])
s8 = gmsh.model.geo.addPlaneSurface([cl9])
s9 = gmsh.model.geo.addPlaneSurface([cl10])

number = 10

# Mesh constraints using transfinite meshing
gmsh.model.geo.mesh.setTransfiniteCurve(l1, number)
gmsh.model.geo.mesh.setTransfiniteCurve(l2, number)
gmsh.model.geo.mesh.setTransfiniteCurve(l3, number)
gmsh.model.geo.mesh.setTransfiniteCurve(l4, number)
gmsh.model.geo.mesh.setTransfiniteCurve(l5, number)
gmsh.model.geo.mesh.setTransfiniteCurve(l6, number)
gmsh.model.geo.mesh.setTransfiniteCurve(l7, number)
gmsh.model.geo.mesh.setTransfiniteCurve(l8, number)
gmsh.model.geo.mesh.setTransfiniteCurve(l9, number)
gmsh.model.geo.mesh.setTransfiniteCurve(l10, number)
gmsh.model.geo.mesh.setTransfiniteCurve(l11, number)
gmsh.model.geo.mesh.setTransfiniteCurve(l12, number)
gmsh.model.geo.mesh.setTransfiniteCurve(l13, number)
gmsh.model.geo.mesh.setTransfiniteCurve(l14, number)
gmsh.model.geo.mesh.setTransfiniteCurve(l15, number)
gmsh.model.geo.mesh.setTransfiniteCurve(l16, number)

gmsh.model.geo.mesh.setTransfiniteCurve(arc1,number)
gmsh.model.geo.mesh.setTransfiniteCurve(arc2,number)
gmsh.model.geo.mesh.setTransfiniteCurve(arc3, number)
gmsh.model.geo.mesh.setTransfiniteCurve(arc4, number)
gmsh.model.geo.mesh.setTransfiniteCurve(arc5, number)
gmsh.model.geo.mesh.setTransfiniteCurve(arc6, number)
gmsh.model.geo.mesh.setTransfiniteCurve(arc7, number)
gmsh.model.geo.mesh.setTransfiniteCurve(arc8, number)


gmsh.model.geo.mesh.setTransfiniteSurface(s2)
gmsh.model.geo.mesh.setTransfiniteSurface(s3)
gmsh.model.geo.mesh.setTransfiniteSurface(s4)
gmsh.model.geo.mesh.setTransfiniteSurface(s5)
gmsh.model.geo.mesh.setTransfiniteSurface(s6)
gmsh.model.geo.mesh.setTransfiniteSurface(s7)
gmsh.model.geo.mesh.setTransfiniteSurface(s8)
gmsh.model.geo.mesh.setTransfiniteSurface(s9)

gmsh.model.geo.mesh.setRecombine(2, s2)
gmsh.model.geo.mesh.setRecombine(2, s3)
gmsh.model.geo.mesh.setRecombine(2, s4)
gmsh.model.geo.mesh.setRecombine(2, s5)
gmsh.model.geo.mesh.setRecombine(2, s6)
gmsh.model.geo.mesh.setRecombine(2, s7)
gmsh.model.geo.mesh.setRecombine(2, s8)
gmsh.model.geo.mesh.setRecombine(2, s9)

# Synchronize and generate mesh
gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)

'''
# Add physical groups for curves and surfaces
gmsh.model.addPhysicalGroup(1, [l9, l10, l11, l12], 11 , name="outer")
gmsh.model.addPhysicalGroup(1, [l1, l2, l3, l4], 12 , name="inner")
gmsh.model.addPhysicalGroup(2, [s1, s2], 13 , name="fluid")
'''
# Save and finalize
gmsh.write("grad_baba_alayina_gider_tek_yeter06.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
# This should be called when you are done using the Gmsh Python API:
gmsh.finalize()
