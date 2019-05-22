import svgwrite
from math import *
from numpy import interp

class Plotdraw:
    def __init__ (self):
        self.paths = []
        self._stroke = [0,0,0]
        self.offsetX = 0
        self.offsetY = 0
        self.mutateExp = "[x,y]"

    def stroke(self, r, g, b):
        self._stroke = [r, g, b]
    
    def line(self, x0, y0, x1, y1):
        x0 += self.offsetX
        y0 += self.offsetY
        x1 += self.offsetX
        y1 += self.offsetY
        self.paths.append([self._stroke,[[x0, y0],[x1, y1]]])

    def linePoints(self, xy0, xy1, points):
        self.beginPath()
        for p in range(0, points):
            x = interp(p, [0, points], [xy0[0], xy1[0]])
            y = interp(p, [0, points], [xy0[1], xy1[1]])
            self.addPoint(x, y)

    def beginPath(self):
        self.paths.append([self._stroke,[]])

    def closePath(self):
        self.paths[len(self.paths)-1][1].append(self.paths[len(self.paths)-1][1][0])


    def addPoint(self, x, y):
        x += self.offsetX
        y += self.offsetY
        xy = self.mutate(x, y)
        x = xy[0]
        y = xy[1]
        self.paths[len(self.paths)-1][1].append([x,y])
    
    def translate(self, x, y):
        self.offsetX += x
        self.offsetY += y

    def setMutate(self, exp):
        self.mutateExp = exp

    def mutate(self, x, y):
        xy = eval(self.mutateExp)
        return xy

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
        
        x = r * sin(radians(0+rot))
        y = r * cos(radians(0+rot))
        self.addPoint(xoff+x, yoff+y)


    def circle(self, xoff, yoff, r):
        
        self.beginPath()

        points = int(pi * (r*2))
        points *= 2 # 1 = 1mm; 2 = 0.5mm; 4 = 0.25mm resolution
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

        for d in range(int(fromAng), int(toAng)):
            x = r * sin(radians(d))
            y = r * cos(radians(d))
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
        self.beginPath()
        for t in range(0, points+1):
            print(t/points)
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
        print(points)
        self.beginPath()
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
        print(points)
        self.beginPath()
        for t in range(0, points+1):
            xy = self.calcBezierQuad(t/points, p0, p1, p2)
            self.addPoint(xy[0], xy[1])

    def bezierQuadPoints(self, p0, p1, p2, points):
        self.beginPath()
        for t in range(0, points+1):
            print(t/points)
            xy = self.calcBezierQuad(t/points, p0, p1, p2)
            self.addPoint(xy[0], xy[1])

    def calcBezierQuad(self, t, p0, p1, p2):
        x = p0[0]*t*t + p1[0]*2*t*(1-t) + p2[0]*(1-t)*(1-t)
        y = p0[1]*t*t + p1[1]*2*t*(1-t) + p2[1]*(1-t)*(1-t)
        return [x,y]

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
                                #join the two paths together and insert it where l1 is
                                ##### try just extending the path — quicker
                                idx = self.paths.index(p1)
                                self.paths.remove(p1)
                                self.paths.remove(p2)
                                self.paths.insert(idx, [p1[0], p1[1] + p2[1][1:]])
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
    
    def distToLine(self, tp,p1,p2):
        return abs((tp[0]-p2[0])*(p2[1]-p1[1]) - (p2[0]-p1[0])*(tp[1]-p2[1])) / sqrt((tp[0]-p2[0])**2 + (tp[1]-p2[1])**2)

    def optimise(self):
        sortedPaths = []
        currPath = self.paths[0]
        while(len(self.paths) > 0):
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


    def HPGL(self, filename, pen, speed):
        scale = 1/0.025
        hpgl  = ""
        hpgl += "IN;VS{};SP{};\n".format(speed, pen)
        for p in self.paths:
            hpgl += "PU{},{};".format(int(p[1][0][0]*scale),int(p[1][0][1]*scale))
            hpgl += "PD"
            for c in p[1][1:-1]:
                hpgl += "{},{},".format(int(c[0]*scale), int(c[1]*scale))
            hpgl += "{},{};\n".format(int(p[1][-1][0]*scale), int(p[1][-1][1]*scale))
        hpgl += "SP0;PU0,0;IN;"

        f = open(filename, "w+")
        f.write(hpgl)
        f.close()

        return hpgl

    def SVG(self, filename, width, height):
        svg = svgwrite.Drawing(filename, profile='full', size=(width, height), debug=False)
        svg.add(svg.rect(insert=(0, 0), size=('100%', '100%'), rx=None, ry=None, fill='rgb(255,255,255)'))


        for p in self.paths:
            path = svg.path(d="M {},{}".format(p[1][0][0],p[1][0][1]))
            for c in p[1][1:]:
                #print(c)
                path.push("L {},{}".format(c[0],c[1]))
            path.stroke(color="rgb({},{},{})".format(p[0][0],p[0][1],p[0][2]))
            path.stroke(width="0.5")
            path.fill(color="none")
            svg.add(path)
        
        svg.save()

        return svg.tostring()

