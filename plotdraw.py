import svgwrite
from math import *
from numpy import interp

class Plotdraw:
    def __init__ (self):
        self.paths = []
        self._stroke = [0,0,0]
        # A2=594,420 - A3=420,297
        self.width = 594
        self.height = 420
        # DPX-2200 (A2) = [-309, -210]
        # DXY-XXXX = [0, 0]
        self.plotOffset = [-309,-210]
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

    def stroke(self, r, g, b):
        self._stroke = [r, g, b]
    
    def line(self, x0, y0, x1, y1):
        self.beginPath()
        self.addPoint(x0, y0)
        self.addPoint(x1, y1)
    
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
        self.paths.append([self._stroke,[]])

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

        for d in range(0,360,int(360/sides)):
            x = r * sin(radians(d+rot))
            y = r * cos(radians(d+rot))
            self.addPoint(xoff+x, yoff+y)
        
        self.closePath()

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

    def bezierPoints(self, p0, p1, p2, p3, points):
        for t in range(0, points+1):
            xy = self.calcBezier(t/points, p0, p1, p2, p3)
            self.addPoint(xy[0], xy[1])

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

    def bezierQuadPoints(self, p0, p1, p2, points):
        for t in range(0, points+1):

            xy = self.calcBezierQuad(t/points, p0, p1, p2)
            self.addPoint(xy[0], xy[1])

    def calcBezierQuad(self, t, p0, p1, p2):
        x = p0[0]*t*t + p1[0]*2*t*(1-t) + p2[0]*(1-t)*(1-t)
        y = p0[1]*t*t + p1[1]*2*t*(1-t) + p2[1]*(1-t)*(1-t)
        return [x,y]

    #returns coordinate of intersection
    def findIntersection(self, x1,y1,x2,y2,x3,y3,x4,y4):
        px = ( (x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4) ) / ( (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4) ) 
        py = ( (x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4) ) / ( (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4) )
        return [px, py]

    def roundValues(self, amount, sigfigs):
        for path in self.paths:
            for p in path[1]:
                for i, v in enumerate(p):
                    p[i] = round((round (p[i] / amount) * amount), sigfigs)

    def concate(self):
        joinCount = 0
        found = True
        #keep looping until no more joins can be found
        while (found == True):
            found = False
            for p1i, p1 in enumerate(self.paths):
                for p2i, p2 in enumerate(self.paths):
                    #don't compare a path with itself
                    if (p1i != p2i):
                        #only compare if they're the same colour
                        if (p1[0] == p2[0]):
                            contig = False
                            #compare the last element of p1 with the first of p2
                            if (p1[1][-1] == p2[1][0]):
                                contig = True
                            #compare the last element of p1 with the last of p2 and reverse if the same
                            if (p1[1][-1] == p2[1][-1]):
                                p2[1].reverse()
                                contig = True
                            #compare the first element of p1 with the last element of p2 and reverse both if same
                            if (p1[1][0] == p2[1][-1]):
                                p1[1].reverse()
                                p2[1].reverse()
                                contig = True
                            #compare the first element of p1 and the first element of p2 and reverse p1 if same
                            if (p1[1][0] == p2[1][0]):
                                p1[1].reverse()
                                contig = True
                            #if contiguous
                            if contig:
                                #join p1 and p2 and remove p2
                                p1.extend(p2[1][1:])
                                self.paths.remove(p2)
                                joinCount += 1
                                found = True
                                break
        return joinCount

    def simplify(self):
        for path in self.paths:
            if path[1][0] != path[1][-1]:
                if len(path[1]) > 2:
                    simplify = True
                    for c in path[1][1:-1]:
                        if self.distToLine(c, path[1][0], path[1][-1]) > 0.5 : simplify = False
                    if simplify : path[1] = [path[1][0],path[1][-1]]

    def outsideBounds(self, xyt, xy0, xy1):
        return (xyt[0] < xy0[0]) or (xyt[0] > xy1[0]) or (xyt[1] < xy0[1]) or (xyt[1] > xy1[1])

    def trimPaths(self, xy0, xy1):
        for i in range(len(self.paths)):
            path = self.paths.pop(0)
            outside = 0
            for xy in path[1]:
                if self.outsideBounds(xy, xy0, xy1):
                    outside = 1
            if outside != 1:
                self.paths.append(path)

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

    def optimise(self):
        print("Optimising...")
        sortedPaths = []
        currPath = self.paths[0]
        while(len(self.paths) > 0):
            print(len(self.paths))
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

    def overplot(self):
        for path in self.paths:
            #reverse the list and remove the first value, and add to the path
            path[1].extend(path[1][::-1][1:])

    def HPGL(self, filename, pen, speed, force):
        scale = 1/0.025
        hpgl  = ""
        hpgl += "IN;VS{};SP{};FS{};\n".format(speed, pen, force)
        for p in self.paths:
            hpgl += "PU{},{};".format(int((p[1][0][0]+self.plotOffset[0])*scale),int((p[1][0][1]+self.plotOffset[1])*scale))
            hpgl += "PD"
            for c in p[1][1:-1]:
                hpgl += "{},{},".format(int((c[0]+self.plotOffset[0])*scale), int((c[1]+self.plotOffset[1])*scale))
            hpgl += "{},{};\n".format(int((p[1][-1][0]+self.plotOffset[0])*scale), int((p[1][-1][1]+self.plotOffset[1])*scale))
        hpgl += "SP0;IN;"

        f = open(filename, "w+")
        f.write(hpgl)
        f.close()

        return hpgl

    def SVG(self, filename):
        svg = svgwrite.Drawing(filename, profile='full', size=(self.width, self.height), debug=False)
        svg.add(svg.rect(insert=(0, 0), size=('100%', '100%'), rx=None, ry=None, fill='rgb(255,255,255)'))


        for p in self.paths:
            path = svg.path(d="M {},{}".format(p[1][0][0],p[1][0][1]))
            for c in p[1][1:]:
                path.push("L {},{}".format(c[0],c[1]))
            path.stroke(color="rgb({},{},{})".format(p[0][0],p[0][1],p[0][2]))
            path.stroke(width="0.5")
            path.fill(color="none")
            svg.add(path)
        
        svg.save()

        return svg.tostring()

