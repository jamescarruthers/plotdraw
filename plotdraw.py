import math, svgwrite

class Plotdraw:
    def __init__ (self):
        self.paths = []
        self._stroke = [0,0,0]

    def stroke(self, r, g, b):
        self._stroke = [r, g, b]
    
    def line(self, x0, y0, x1, y1):
        self.paths.append([self._stroke,[[x0, y0],[x1, y1]]])

    def beginPath(self):
        self.paths.append([self._stroke,[]])
        
    def addPoint(self, x, y):
        self.paths[len(self.paths)-1][1].append([x,y])
        
    def poly(self, xoff, yoff, r, sides, rot):
        
        self.beginPath()

        for d in range(0,360,int(360/sides)):
            x = r * math.sin(math.radians(d+rot))
            y = r * math.cos(math.radians(d+rot))
            self.addPoint(xoff+x, yoff+y)
        
        x = r * math.sin(math.radians(0+rot))
        y = r * math.cos(math.radians(0+rot))
        self.addPoint(xoff+x, yoff+y)

    def bezier(self, p0, p1, p2, p3, points):
        self.beginPath()
        for t in range(0, points+1):
            print(t/points)
            xy = self.calcBezier(t/points, p0, p1, p2, p3)
            self.addPoint(xy[0], xy[1])

    def calcBezier(self, t, p0, p1, p2, p3):
        x = (1-t)*(1-t)*(1-t)*p0[0] + 3*(1-t)*(1-t)*t*p1[0] + 3*(1-t)*t*t*p2[0] + t*t*t*p3[0]
        y = (1-t)*(1-t)*(1-t)*p0[1] + 3*(1-t)*(1-t)*t*p1[1] + 3*(1-t)*t*t*p2[1] + t*t*t*p3[1]
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
        return abs((tp[0]-p2[0])*(p2[1]-p1[1]) - (p2[0]-p1[0])*(tp[1]-p2[1])) / math.sqrt((tp[0]-p2[0])**2 + (tp[1]-p2[1])**2)

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
        return sorted([[abs(math.sqrt(pow((x[0]-coord[0]),2) + pow((x[1]-coord[1]),2))),i] for i, x in enumerate(coords)])[0]


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

