
"""
"LBLABEL COMMAND";CHR$(3) - lower left origin - terminate with character ASCII 03 ETX
"SI0.3,0.4;" - CHAR SIZE IN CM
"SL0.4;" - SLANT CHARACTER
"""

import svgwrite
from math import *
from numpy import interp
import noise
import time

class Plotdraw:
    def __init__ (self):
        # path format
        # [[[pen,speed,force],[fill,hatch,angle]],[[x,y],[x,y],[x,y]]]
        self.paths = []
        self.penOptions = [1, 24, 4] # [pen, speed, force]
        # lineoption - [0]- 0 = no line (polygon only), 1 = line, 2 = pattern (uses [1] and [2] as options)
        self.lineOptions = [1, 1, 4] # [line, pattern, dot/dash pitch]
        self.fillOptions = [0, 2, 45] # [pattern, spacing, angle]
        # A2=594,420 - A3=420,297
        self.width = 594
        self.height = 420
        # DPX-2200 (A2) = [-309, -210]
        # DXY-XXXX = [0, 0]
        self.plotOffset = [-309,-210]
        self.inputWindow = []
        self.offsetX = 0
        self.offsetY = 0
        self.mutateX = "x"
        self.mutateY = "y"
        self.d = 50
        self.perspective = self.width * 0.8
        self.isoAngle = 30
        self.zAngle = 0
        self.rotateCentre = [0,0,0]
        self.rotateInc = 0
        self.rotateAzi = 0
        # [[x,y],[w,h],title]
        self.title = []

    def pen(self, pen):
        self.penOptions[0] = pen

    def penSpeed(self, speed):
        self.penOptions[1] = speed

    def penForce(self, force):
        self.penOptions[2] = force

    # Fill types
    # 0: no fill
    # 1: cross directional fill (pen thickness)
    # 2: one directional fill (pen thickness)
    # 3: hatching (spacing)
    # 4: cross hatching (spacing)
    # 5: cross directional hatching (user defined fill)
    # 6: one directional hatching (user defined fill)
    def fill(self, pattern):
        self.fillOptions[0] = pattern

    # spacing between hatched lines
    def fillSpacing(self, spacing):
        self.fillOptions[1] = spacing

    # angle of hatching in degrees
    def fillAngle(self, angle):
        self.fillOptions[2] = angle

    # 0 = no line, 1 = solid line, 2 = pattern line    
    def lineType(self, line):
        self.lineOptions[0] = line

    # 0 = dot on coordinates 1 . . . 2 - - - 3 — — — 4 —.—.—. 5 —-—-—- 6 —--—--—-- 
    def linePattern(self, pattern):
        self.lineOptions[1] = pattern

    # percentage of distance between P1 and P2
    def linePitch(self, pitch):
        self.lineOptions[2] = pitch

    def line(self, xy0, xy1):
        self.beginPath()
        self.addPoint(xy0[0], xy0[1])
        self.addPoint(xy1[0], xy1[1])
    
    def line3d(self, x0, y0, z0, x1, y1, z1):
        self.beginPath()
        self.addPoint3d(x0, y0, z0)
        self.addPoint3d(x1, y1, z1)

    def linePoints(self, xy0, xy1, points):
        self.beginPath()
        for p in range(0, points):
            x = interp(p, [0, points], [xy0[0], xy1[0]])
            y = interp(p, [0, points], [xy0[1], xy1[1]])
            self.addPoint(x, y)

    def linePoints3d(self, xyz0, xyz1, points):
        self.beginPath()
        for p in range(0, points):
            x = interp(p, [0, points], [xyz0[0], xyz1[0]])
            y = interp(p, [0, points], [xyz0[1], xyz1[1]])
            z = interp(p, [0, points], [xyz0[2], xyz1[2]])
            self.addPoint3d(x, y, z)

    def beginPath(self):
        self.paths.append([[self.penOptions[:], self.lineOptions[:], self.fillOptions[:]],[]])

    def closePath(self):
        self.paths[len(self.paths)-1][1].append(self.paths[len(self.paths)-1][1][0])

    def addPoint(self, x, y):
        x += self.offsetX
        y += self.offsetY
        x = eval(self.mutateX)
        y = eval(self.mutateY)
        self.paths[len(self.paths)-1][1].append([x,y])
    
    def addPoint3d(self, x, y, z):
        x,y,z = self.rotate3d(x,y,z)
        scalePerspective = self.perspective / (self.perspective + z)
        x = x * scalePerspective
        y = y * scalePerspective
        self.addPoint(x, y)

    def addPointIso(self, x, y, z):
        x,y,z = self.rotate3d(x,y,z)
        xp = (x * cos(radians(-self.isoAngle))) - (y * cos(radians(-180+self.isoAngle)))
        yp = z + (x * sin(radians(-self.isoAngle))) - (y * sin(radians(-180+self.isoAngle)))
        self.addPoint(xp, yp)

    def rotate3d(self, x, y, z):
        rad, inc, azi = self.cartToSph(x-self.rotateCentre[0],y-self.rotateCentre[1],z-self.rotateCentre[2])
        
        inc += self.rotateInc
        azi += self.rotateAzi

        x,y,z = self.sphToCart(rad, inc, azi)

        #xp = x*cos(q) - y*sin(q)
        #yp = x*sin(q) + y*cos(q)
        #zp = z

        return x+self.rotateCentre[0],y+self.rotateCentre[1],z+self.rotateCentre[2]

    def sphToCart(self, rho, phi, theta):
        x = rho * sin(radians(phi)) * cos(radians(theta))
        y = rho * sin(radians(phi)) * sin(radians(theta))
        z = rho * cos(radians(phi))
        return x, y, z

    def cartToSph2(self, x, y, z):
        rad = sqrt(x*x + y*y + z*z)
        if rad != 0 : inc = degrees(acos(z/rad))
        if rad == 0 : inc = 0
        azi = degrees(atan2(y,x))
        return rad, inc, azi

    def cartToSph(self, x, y, z):
        rho = sqrt(x**2 + y**2 + z**2)
        theta = atan2(y, x)
        if rho != 0.0:
            phi = acos( z / rho )
        else:
            #phi = pi / 2 * (z)
            phi = 0
        return rho, degrees(phi), degrees(theta)

    #run every path through the mutate expression
    def mutate(self):
        for path in self.paths:
            for xy in path[1]:
                x = xy[0]
                y = xy[1]
                x = eval(self.mutateX)
                y = eval(self.mutateY)
                xy[0] = x
                xy[1] = y

    def translate(self, x, y):
        self.offsetX += x
        self.offsetY += y

    def straight(self, x, y):
        return [x,y]

    def offset(self, x, y):
        self.offsetX = x
        self.offsetY = y

    def poly(self, xoff, yoff, r, sides, rot):
        
        self.beginPath()

        for side in range(0, sides):
            d = side * (360/sides)
            x = r * sin(radians(d + rot))
            y = r * cos(radians(d + rot))
            self.addPoint(xoff + x, yoff + y)
        
        self.closePath()

    def polyCircum(self, xoff, yoff, r, sides, rot):

        #calculate radius
        r = r / cos(radians(360/(sides*2)))

        self.poly(xoff, yoff, r, sides, rot)
        

    def polyPinch(self, xoff, yoff, r, sides, rot, pinch):
        
        self.beginPath()

        deg = 360 / sides

        for s in range(0, sides):

            fromDeg = (s * deg) + rot
            toDeg = ((s+1) * deg) + rot

            if (isinstance(pinch, list)):
                rPinch = r * pinch[s % len(pinch)]
            else:
                rPinch = r * pinch
            
            p0 = [r * sin(radians(fromDeg)),r * cos(radians(fromDeg))]
            p1 = [(rPinch) * sin(radians(fromDeg)),(rPinch) * cos(radians(fromDeg))]
            p2 = [(rPinch) * sin(radians(toDeg)),(rPinch) * cos(radians(toDeg))]
            p3 = [r * sin(radians(toDeg)),r * cos(radians(toDeg))]

            self.bezier(p0, p1, p2, p3)
            #self.line(p0[0],p0[1],p3[0],p3[1])

    def circle2p(self, xy0, xy1):

        midpoint = [(xy0[0]+xy1[0])/2,(xy1[1]+xy1[1])/2]
        radius = sqrt(pow((xy1[0]-xy0[0]), 2) + pow((xy1[1]-xy0[1]), 2)) / 2
        self.circle(midpoint[0], midpoint[1], radius)


    def circle(self, xoff, yoff, r):

        r = abs(r)
        
        self.beginPath()

        points = int(pi * (r*2))
        points *= 2 # 1 = 1mm; 2 = 0.5mm; 4 = 0.25mm resolution
        if points < 4 : points = 4
        pointsInc = 360 / points

        for d in range(0, points):
            x = r * sin(radians(d*pointsInc))
            y = r * cos(radians(d*pointsInc))
            self.addPoint(xoff+x, yoff+y)
        
        x = r * sin(radians(0))
        y = r * cos(radians(0))
        self.addPoint(xoff+x, yoff+y)

    def arc(self, xoff, yoff, r, fromAng, toAng):

        #self.beginPath()

        angle = abs((fromAng-toAng)%360)
        distance = int(2 * pi * r * (angle/360.0))
        #distance is in mm, so multiply to add resolution — 2x = 0.5mm, 4x = 0.25mm

        if distance == 0: distance = 1

        for d in range(0, distance):
            deg = interp(d, [0, distance], [fromAng, toAng])
            x = r * sin(radians(deg))
            y = r * cos(radians(deg))
            self.addPoint(xoff+x, yoff+y)


    def spiral(self, xoff, yoff, fromAng, toAng, a, b):

        self.beginPath()
        for ang in range(fromAng, toAng):
            r = a + (b*radians(ang))
            x = r * cos(radians(ang))
            y = r * sin(radians(ang))
            self.addPoint(xoff+x, yoff+y)


    def polyRound(self, xoff, yoff, r, sides, rot, rnd):

        r = r - rnd
        ang = 360/sides
        halfang = ang/2

        self.beginPath()

        for i in range(0, sides):
            angle = i*ang+rot
            x = r * sin(radians(angle))
            y = r * cos(radians(angle))
            self.arc(x,y,rnd,angle-halfang,angle+halfang)

        self.closePath()


    def polyRound2(self, xoff, yoff, r, sides, rot, rnd):
        
        self.beginPath()

        deg = 360 / sides

        for s in range(0, sides):

            if (isinstance(rot, list)):
                thisRot = rot[s % len(rot)]
            else:
                thisRot = rot

            fromDeg = (s * deg) + thisRot
            toDeg = ((s+1) * deg) + thisRot

            if (isinstance(r, list)):
                thisR = r[s % len(r)]
            else:
                thisR = r

            if (isinstance(rnd, list)):
                thisRnd = thisR * rnd[s % len(rnd)]
            else:
                thisRnd = thisR * rnd

            p0 = [thisR * sin(radians(fromDeg)), thisR * cos(radians(fromDeg))]
            p1 = [p0[0]+(thisRnd * sin(radians(fromDeg+deg))), p0[1]+(thisRnd * cos(radians(fromDeg+deg)))]
            p3 = [thisR * sin(radians(toDeg)),thisR * cos(radians(toDeg))]
            p2 = [p3[0]+(thisRnd * sin(radians(toDeg-deg))), p3[1]+(thisRnd * cos(radians(toDeg-deg)))]

            self.bezier(p0, p1, p2, p3)
            #self.line(p0[0],p0[1],p3[0],p3[1])

        self.closePath()


    def ellipse(self, xoff, yoff, rw, rh):
        
        self.beginPath()

        r = max([rw,rh])
        points = int(pi * (r*2))
        points *= 2 # 1 = 1mm; 2 = 0.5mm; 4 = 0.25mm resolution
        pointsInc = 360 / points

        for d in range(0,points):
            x = rw * sin(radians(d*pointsInc))
            y = rh * cos(radians(d*pointsInc))
            self.addPoint(xoff+x, yoff+y)
        
        x = rw * sin(radians(0))
        y = rh * cos(radians(0))
        self.addPoint(xoff+x, yoff+y)

    #draw a bezier with a set number of points
    def bezierPoints(self, p0, p1, p2, p3, points):
        for t in range(0, points+1):
            xy = self.calcBezier(t/points, p0, p1, p2, p3)
            self.addPoint(xy[0], xy[1])

    #draw a bezier
    def bezier(self, p0, p1, p2, p3):
        distance = 0
        pointsTest = 8
        for t in range(0, pointsTest-1):
            xy0 = self.calcBezier(t/pointsTest, p0, p1, p2, p3)
            xy1 = self.calcBezier((t+1)/pointsTest, p0, p1, p2, p3)
            distance += abs(sqrt(((xy1[0]-xy0[0])**2)+((xy1[1]-xy0[1])**2)))
        
        points = int(distance*1) #multiplier to add or decrease resolution

        for t in range(0, points+1):
            xy = self.calcBezier(t/points, p0, p1, p2, p3)
            self.addPoint(xy[0], xy[1])


    def calcBezier(self, t, p0, p1, p2, p3):
        x = (1-t)*(1-t)*(1-t)*p0[0] + 3*(1-t)*(1-t)*t*p1[0] + 3*(1-t)*t*t*p2[0] + t*t*t*p3[0]
        y = (1-t)*(1-t)*(1-t)*p0[1] + 3*(1-t)*(1-t)*t*p1[1] + 3*(1-t)*t*t*p2[1] + t*t*t*p3[1]
        return [x,y]

    #draw a quad bezier
    def bezierQuad(self, p0, p1, p2):
        distance = 0
        pointsTest = 8
        for t in range(0, pointsTest-1):
            xy0 = self.calcBezierQuad(t/pointsTest, p0, p1, p2)
            xy1 = self.calcBezierQuad((t+1)/pointsTest, p0, p1, p2)
            distance += abs(sqrt(((xy1[0]-xy0[0])**2)+((xy1[1]-xy0[1])**2)))
        
        points = int(distance*1) #multiplier to add or decrease resolution

        for t in range(0, points+1):
            xy = self.calcBezierQuad(t/points, p0, p1, p2)
            self.addPoint(xy[0], xy[1])

    #draw a quad bezier with a set number of points
    def bezierQuadPoints(self, p0, p1, p2, points):
        for t in range(0, points+1):
            xy = self.calcBezierQuad(t/points, p0, p1, p2)
            self.addPoint(xy[0], xy[1])

    def calcBezierQuad(self, t, p0, p1, p2):
        x = p0[0]*t*t + p1[0]*2*t*(1-t) + p2[0]*(1-t)*(1-t)
        y = p0[1]*t*t + p1[1]*2*t*(1-t) + p2[1]*(1-t)*(1-t)
        return [x,y]

    #returns coordinate of intersection
    def lineIntersection(self, p0, p1, p2, p3):
        Bx_Ax = p1[0] - p0[0]
        By_Ay = p1[1] - p0[1]
        Dx_Cx = p3[0] - p2[0]
        Dy_Cy = p3[1] - p2[1]
        determinant = (-Dx_Cx * By_Ay + Bx_Ax * Dy_Cy) 
        if abs(determinant) < 1e-20: 
            return None 
        s = (-By_Ay * (p0[0] - p2[0]) + Bx_Ax * (p0[1] - p2[1])) / determinant 
        t = ( Dx_Cx * (p0[1] - p2[1]) - Dy_Cy * (p0[0] - p2[0])) / determinant 
        if s >= 0 and s <= 1 and t >= 0 and t <= 1: 
            return [p0[0] + (t * Bx_Ax), p0[1] + (t * By_Ay)]
        return None

    def sortXY(self, val):
        return val[0]+val[1]

    def hatch(self, angle, spacing):

        paths = [self.paths[-1], self.paths[-2]] 

        #minxy = [min([xy[0] for xy in path[1]]), min([xy[1] for xy in path[1]])]
        #maxxy = [max([xy[0] for xy in path[1]]), max([xy[1] for xy in path[1]])]
        
        # should really be furthest ll point and furthest ur point in path

        minxy = min(min([path[1] for path in paths]))
        maxxy = max(max([path[1] for path in paths]))

        #minxy = [0,0]
        #maxxy = [self.width, self.height]
        
        midx = (minxy[0] + maxxy[0]) / 2
        midy = (minxy[1] + maxxy[1]) / 2
        diameter = abs(self.dist(minxy, maxxy))
        
        radius = (diameter/2) / cos(radians(360/8))
        
        lx = midx + (radius * cos(radians(45+angle)))
        ly = midy + (radius * sin(radians(45+angle)))

        rx = midx + (radius * cos(radians(45+90+angle)))
        ry = midy + (radius * sin(radians(45+90+angle)))

        distance = self.dist([lx, ly],[rx, ry])

        xoff = spacing * cos(radians(90+angle))
        yoff = spacing * sin(radians(90+angle))

        number = ceil(distance / spacing)
        
        hatched = []
        
        for i in range(0, number):
            ixoff = xoff * i
            iyoff = yoff * i
            if i % 2 == 0:
                #self.line([lx-ixoff, ly-iyoff],[rx-ixoff, ry-iyoff])
                hatched.append([[lx-ixoff, ly-iyoff],[rx-ixoff, ry-iyoff]])
            else:
                #self.line([lx-ixoff, ly-iyoff],[rx-ixoff, ry-iyoff])
                hatched.append([[rx-ixoff, ry-iyoff],[lx-ixoff, ly-iyoff]])

        for line in hatched:
            for h0, h1 in zip(line, line[1:]):
                newLine = []
                for path in paths:
                    for p0, p1 in zip(path[1], path[1][1:]):
                        inter = self.lineIntersection(h0, h1, p0, p1)
                        if inter != None:
                            newLine.append(inter)
                    if len(newLine) > 1 and len(newLine) % 2 == 0:
                        print("unsorted ({}): ".format(len(newLine)), end="")
                        print(newLine)
                        #newLine = sorted(newLine, key=lambda k: [k[1], k[0]])
                        newLine.sort()
                        print("sorted: ", end="")
                        print(newLine)
                        for i in range(0, len(newLine), 2):
                            self.beginPath()
                            self.paths[len(self.paths)-1][1].append([newLine[i][0], newLine[i][1]])
                            self.paths[len(self.paths)-1][1].append([newLine[i+1][0], newLine[i+1][1]])
        

    def roundValues(self, amount, sigfigs):
        for path in self.paths:
            for p in path[1]:
                for i, v in enumerate(p):
                    p[i] = round((round (p[i] / amount) * amount), sigfigs)

    #connect conjoining paths
    def concate(self):
        startTime = time.time()
        startCount = len(self.paths)
        print("Concating {} paths".format(len(self.paths)))
        found = True
        while (found == True):
            found = False
            for p1 in self.paths:
                for p2 in self.paths:
                    if p1 != p2:
                        join = False
                        # end of p1 == start of p2
                        if p1[1][-1] == p2[1][0]:
                            join = True
                        # end of p1 == end of p2
                        elif p1[1][-1] == p2[1][-1]:
                            p2[1].reverse()
                            join = True
                        # start of p1 == start of p2
                        elif p1[1][0] == p2[1][0]:
                            p1[1].reverse()
                            join = True
                        # start of p1 == end of p2
                        elif p1[1][0] == p2[1][-1]:
                            p1[1].reverse()
                            p2[1].reverse()
                        if join:
                            p1[1].extend(p2[1][1:])
                            self.paths.remove(p2)
                            found = True
            print("Reducing to {} paths".format(len(self.paths)), end = '\r')

        print()
        print("Joined {} paths".format(startCount - len(self.paths)))
        print(time.time() - startTime)

    #remove points from paths that are too close to each other
    def simplify(self):
        for path in self.paths:
            if path[1][0] != path[1][-1]:
                if len(path[1]) > 2:
                    simplify = True
                    for c in path[1][1:-1]:
                        if self.distToLine(c, path[1][0], path[1][-1]) > 0.5 : simplify = False
                    if simplify : path[1] = [path[1][0],path[1][-1]]
    
    #used in trim functions
    def outsideBounds(self, xyt, xy0, xy1):
        return (xyt[0] < xy0[0]) or (xyt[0] > xy1[0]) or (xyt[1] < xy0[1]) or (xyt[1] > xy1[1])

    #remove paths outside a bounding box
    def trimPaths(self, xy0, xy1):
        for i in range(len(self.paths)):
            path = self.paths.pop(0)
            outside = 0
            for xy in path[1]:
                if self.outsideBounds(xy, xy0, xy1):
                    outside = 1
            if outside != 1:
                self.paths.append(path)

    #remove points outside a bounding box
    def trimPoints(self, xy0, xy1):
        for i in range(len(self.paths)):
            path = self.paths.pop(0)
            inside = 0
            for xy in path[1]:
                if self.outsideBounds(xy, xy0, xy1):
                    inside = 0
                else:
                    if inside == 0:
                        self.beginPath()
                        inside = 1
                    self.addPoint(xy[0], xy[1])

    def distToLine(self, tp,p1,p2):
        return abs((tp[0]-p2[0])*(p2[1]-p1[1]) - (p2[0]-p1[0])*(tp[1]-p2[1])) / sqrt((tp[0]-p2[0])**2 + (tp[1]-p2[1])**2)

    def dist(self, p0, p1):
        return sqrt(pow(p1[0]-p0[0],2)+pow(p1[1]-p0[1],2))

    #optimise path order to reduce pen up time
    def optimise(self):
        print("Optimising...")
        sortedPaths = []
        currPath = self.paths[0]
        while(len(self.paths) > 0):
            #print("{}.. ".format(len(self.paths)))
            if (len(self.paths) == 1):
                sortedPaths.append(currPath)
                self.paths.remove(currPath)
                break
            sortedPaths.append(currPath)
            lastCoords = currPath[1][-1]
            self.paths.remove(currPath)
            startResult = self.findClosest(lastCoords, [y[0] for y in [x[1] for x in self.paths]])
            endResult = self.findClosest(lastCoords,  [y[-1] for y in [x[1] for x in self.paths]])
            if (startResult[0] < endResult[0]):
                closest = startResult[1]
            else:
                closest = endResult[1]
                #reverse path
                self.paths[closest][1].reverse()
            currPath = self.paths[closest]
        self.paths = sortedPaths

    def findClosest(self, coord, coords):
        return sorted([[abs(sqrt(pow((x[0]-coord[0]),2) + pow((x[1]-coord[1]),2))),i] for i, x in enumerate(coords)])[0]

    def countPoints(self):
        count = 0
        for path in self.paths:
            count += len(path[1])
        return count

    #append a reversed version to every path - 1-2-3-4 -> 1-2-3-4-3-2-1
    def overplot(self):
        for path in self.paths:
            #reverse the list and remove the first value, and add to the path
            path[1].extend(path[1][::-1][1:])

    def centre(self):
        minx = []
        miny = []
        maxx = []
        maxy = []
        
        for path in self.paths:
            minx.append(min([xy[0] for xy in path[1]]))
            miny.append(min([xy[1] for xy in path[1]]))
            maxx.append(max([xy[0] for xy in path[1]]))
            maxy.append(max([xy[1] for xy in path[1]]))

        minx = min(minx)
        miny = min(miny)
        maxx = max(maxx)
        maxy = max(maxy)

        width = maxx - minx
        height = maxy - miny



        

    #output a HPGL file
    def HPGL(self, filename):
        
        print("Generating HPGL...")
        
        scale = 1/0.025
        hpgl  = ""
        hpgl += "IN;DF;\n"
        
        lastPen = []
        lastLine = []
                
        if self.inputWindow != []:
            hpgl += "IW{},{},{},{};\n".format(int((self.inputWindow[0][0]+self.plotOffset[0])*scale), int((self.inputWindow[0][1]+self.plotOffset[1])*scale), int((self.inputWindow[1][0]+self.plotOffset[0])*scale), int((self.inputWindow[1][1]+self.plotOffset[1])*scale))
        for p in self.paths:
            # if not a filled polygon
            if (p[0][2][0] == 0):
                if p[0][0] != lastPen:
                    lastPen = p[0][0]
                    hpgl += "SP{};VS{};FS{};\n".format(p[0][0][0], p[0][0][1], p[0][0][2])
                if p[0][1] != lastLine:
                    lastLine = p[0][1]
                    if (p[0][1][0] == 1):
                        hpgl += "LT;\n"
                    if (p[0][1][0] == 2):
                        hpgl += "LT{},{};\n".format(p[0][1][1],p[0][1][2])
                if (p[0][1][0] > 0):
                    hpgl += "PU{},{};".format(int((p[1][0][0]+self.plotOffset[0])*scale),int((p[1][0][1]+self.plotOffset[1])*scale))
                    hpgl += "PD"
                    for c in p[1][1:-1]:
                        hpgl += "{},{},".format(int((c[0]+self.plotOffset[0])*scale), int((c[1]+self.plotOffset[1])*scale))
                    hpgl += "{},{};\n".format(int((p[1][-1][0]+self.plotOffset[0])*scale), int((p[1][-1][1]+self.plotOffset[1])*scale))
            # if a filled polygon
            if (p[0][2][0] > 0):
                if p[0][0] != lastPen:
                    hpgl += "SP{};VS{};FS{};\n".format(p[0][0][0], p[0][0][1], p[0][0][2])
                    lastPen = p[0][0]
                if p[0][1] != lastLine:
                    lastLine = p[0][1]
                    if (p[0][1][0] == 1):
                        hpgl += "LT;\n"
                    if (p[0][1][0] == 2):
                        hpgl += "LT{},{};\n".format(p[0][1][1],p[0][1][2])
                hpgl += "PU{},{};".format(int((p[1][0][0]+self.plotOffset[0])*scale),int((p[1][0][1]+self.plotOffset[1])*scale))
                hpgl += "PM0;PD"
                for c in p[1][1:-1]:
                    hpgl += "{},{},".format(int((c[0]+self.plotOffset[0])*scale), int((c[1]+self.plotOffset[1])*scale))
                hpgl += "{},{};".format(int((p[1][-1][0]+self.plotOffset[0])*scale), int((p[1][-1][1]+self.plotOffset[1])*scale))
                hpgl += "PM2;\n"    
                hpgl += "FT{},{},{};FP;\n".format(p[0][2][0],int(p[0][2][1]*scale),p[0][2][2])
                if (p[0][1][0] > 0):
                    hpgl += "EP;\n"

        hpgl += "SP0;IN;"

        f = open(filename, "w+")
        f.write(hpgl)
        f.close()

        return hpgl

    #output an SVG file
    def SVG(self, filename="output.svg"):
        
        print("Generating SVG...")
        
        svg = svgwrite.Drawing(filename, profile='full', size=(self.width, self.height), debug=False)
        svg.add(svg.rect(insert=(0, 0), size=('100%', '100%'), rx=None, ry=None, fill='rgb(255,255,255)'))

        for p in self.paths:
            path = svg.path(d="M {},{}".format(p[1][0][0],self.height-p[1][0][1]))
            for c in p[1][1:]:
                path.push("L {},{}".format(c[0],self.height-c[1]))
            #path.stroke(color="rgb({},{},{})".format(p[0][0],p[0][1],p[0][2]))
            path.stroke(color="rgb({},{},{})".format(0,0,0))
            path.stroke(width="0.5")
            if (p[0][2][0] == 0):
                path.fill(color="none")
            else:
                path.fill(color="rgb({},{},{})".format(0,0,0))
            svg.add(path)
        
        if self.inputWindow != []:
            window = svg.path(d="M {},{}".format(self.inputWindow[0][0], self.height-self.inputWindow[0][1]))
            window.push("L {},{}".format(self.inputWindow[0][0], self.height-self.inputWindow[1][1]))
            window.push("L {},{}".format(self.inputWindow[1][0], self.height-self.inputWindow[1][1]))
            window.push("L {},{}".format(self.inputWindow[1][0], self.height-self.inputWindow[0][1]))
            window.push("L {},{}".format(self.inputWindow[0][0], self.height-self.inputWindow[0][1]))
            window.stroke(color="rgb({},{},{})".format(0,255,0))
            window.stroke(width="0.2")
            window.fill(color="none")
            svg.add(window)
        svg.save()

        return svg.tostring()

