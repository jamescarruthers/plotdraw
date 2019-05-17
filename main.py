import plotdraw

dwg = plotdraw.Plotdraw()

#define stroke rgb colour.
dwg.stroke(128, 128, 128)

#add two lines. lines are just paths with two points.
print("add lines")
dwg.line(0,0,100,100)
dwg.line(100,100,200,200)
dwg.line(120,120,250,250)
dwg.line(210,210,300,300)
print(dwg.paths)

#add a path by adding new vertex
print("add a path")
dwg.beginPath()
dwg.addPoint(20,30)
dwg.addPoint(40,50)
dwg.addPoint(80,100)
print(dwg.paths)

#add a bezier curve
print("add a path")
dwg.bezier([0,0],[0,100],[400,100],[400,0], 20)
print(dwg.paths)

#5 sided polygon with radius of 40 at 500,400 rotated 0 degrees
print("add a polygon")
dwg.poly(500,400, 40, 5, 0)
print(dwg.paths)

#concate the paths - connect all paths with the same start and end points (iteratively)
print("concate paths")
dwg.concate()
print(dwg.paths)

#simplify the paths - turn paths in to lines if paths points are all on the same line.
print("simplified")
dwg.simplify()
print(dwg.paths)

#work out the most efficient way to physically draw the paths by ordering them so the pen travels the shortest distance in the up position between paths
print("optimised for drawing")
dwg.optimise()
print(dwg.paths)

#this is more for readability and/or might be worth doing before concating and simplifying to connect paths that are close
print("round values to plotter resolution increments")
dwg.roundValues(0.025, 3)
print(dwg.paths)

#output HPGL using a pen 0 and speed of 6mm/s
print("output HPGL")
print(dwg.HPGL("output.hpgl", 0, 6))

#output preview SVG with a width of 420 and height of 300
print("output SVG")
print(dwg.SVG("output.svg", 420, 300))
