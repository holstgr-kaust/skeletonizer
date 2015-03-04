import re


#
# Node class
#

class Node(object):
    """Graph node point class in 3D space """

    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def list(self):
        """ Returns a list of XYZ values"""
        return [self.x, self.y, self.z]

    def position(self):
        """ Returns a tuple of XYZ values"""
        return (self.x, self.y, self.z)

#
# 3D point class
#

class Point3D(Node):
    """3D Point class with public x,y,z attributes and a diameter"""

    def __init__(self, x=0.0, y=0.0, z=0.0, d=0.0):
        self.x = x
        self.y = y
        self.z = z
        self.diameter = d

    def list(self):
        """ Returns a list of XYZD values"""
        return [self.x, self.y, self.z, self.diameter]

    def position(self):
        """ Returns a tuple of XYZ values"""
        return (self.x, self.y, self.z)

    def set_diameter(self, dia):
        """
        Set diameter value
        :param dia: Float diameter value
        """
        self.diameter = dia


#
# Segment class
#

class Segment(object):
    """Array of Point objects together with start and end node names

    Length of a segment is the count of points. End points of a
    segment are at the same location as the end nodes.

    """

    def __init__(self, start, end):
        self.start = start
        self.end = end

        self.pointcount = None
        self.points = []

    def __len__(self):
        """ Count of points"""
        return(self.pointcount)

#
# Skeleton class
#

class Skeleton(object):
    """Top storage object for a skeleton that knows about
    nodes, segments, and their locations """

    def __init__(self):
        self.nodes = {}
        self.segments = []

    def add_node(self, name, node):
        """Add one Node object to a dictionary

        The name is the key to the dictionary of nodes.
        """
        self.nodes.setdefault(name, node)

    def add_segment(self, segment):
        """ Add one Segment object to an array"""
        self.segments.append(segment)

    def add_points(self, points):
        """Add an array of Point objects

        The skeleton needs to be populated with its segments before
        calling this method. Segments need to have the point count
        (Segment.pointcount) for this method to pass the points to
        their correct segments.
        """
        offset = 0
        for segment in self.segments :
            segment.points = points[offset:offset+segment.pointcount]
            offset += segment.pointcount

    def info(self):
        """Print out the count of Node, Segment and Points objects"""
        c = 0
        for s in self.segments:
             c+= len(s.points)
        return "Nodes    : %5i\nSegments : %5i\nPoints   : %5i" % (len(self.nodes), len(self.segments), c)

#
# AmirameshReader class
#

class AmirameshReader(object):
    """ Read from a filehandle, parse, return a Skeleton object"""

    def parse(self, f):

        skel = Skeleton()       # storage object
        points = []             # list of points
        counter = 0             # section counter
        linecounter = 0         # within sections

        for line in f:
            # trim white space, including \r,\n
            line = line.strip()

            # ignore empty lines
            if not line:
                continue

            # skip intro
            header = line.startswith("@")
            if counter == 0 and not header:
                continue

            # header
            if header:
                counter+= 1
                linecounter = 0
                continue

            if counter == 1:            # nodes
                match = re.search('([\d\.e\+\-]+) ([\d\.e\+\-]+) ([\d\.e\+\-]+)', line)
                x,y,z = match.groups()
                x = float(x)
                y = float(y)
                z = float(z)
                n = Node(x,y,z)
                skel.add_node(linecounter,n)
                linecounter += 1

            elif counter == 2:          # segments to nodes
                match = re.search('(\d+) (\d+)', line)
                start,end = match.groups()
                seg = Segment(int(start), int(end))
                skel.add_segment(seg)

            elif counter == 3:          # point count within segment
                match = re.search('(\d+)', line)
                count = match.groups()
                skel.segments[linecounter].pointcount = int(count[0])
                linecounter += 1

            elif counter == 4:          # point coordinates within a segment
                match = re.search('([\d\.e\+\-]+) ([\d\.e\+\-]+) ([\d\.e\+\-]+)', line)
                x,y,z = match.groups()
                x = float(x)
                y = float(y)
                z = float(z)
                p = Point3D(x,y,z)
                points.append(p)
                #linecounter += 1

            elif counter == 5:           # diameter
                # empty values replaced by 0
                if line == "nan":
                    line = "0.0"

                match = re.search('([\d\.e\+\-]+)', line)
                dia = match.groups()[0]
                dia = float(dia)
                points[linecounter].set_diameter(dia)
                linecounter += 1

        # add points in the end for efficiency
        skel.add_points(points)
        return skel
