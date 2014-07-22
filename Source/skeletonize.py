#!/usr/bin/env python

"""
This program creates a BBPSDK HDF5 morphology file from a skeletonization representation in an Amiramesh text file.
"""
import os
import re
import sys
import math
import getopt
import copy
import json
import logging
from collections import defaultdict

from bbp_import_module import *
from amiramesh import *

# TODO: move helper math functions into its own simple_math module, or replace with numpy
def vlogger(func):
    def inner(*args, **kwargs):
        print("Called: %s(%s, %s)" % (func.__name__, args, kwargs))
        result = func(*args, **kwargs)
        print(" result: %s, len:%s\n" % (str(result), str(vlength(result))))
        return result
    return inner

def square(x):
    return x * x

def distance_squared(v1, v2):
    return sum(map(lambda x, y: square(x - y), v1, v2))

def distance_squared(v1, v2):
    return sum(map(lambda x, y: square(x - y), v1, v2))

def distance(v1, v2):
    return math.sqrt(distance_squared(v1, v2))


def vlength(vect):
    return math.sqrt(sum(map(lambda v: square(v), vect)))

def vmuls3(v, x):
    return (v[0]*x, v[1]*x, v[2]*x)

def vdivs3(v, x):
    return (v[0]/x, v[1]/x, v[2]/x)

def vadds3(v, x):
    return (v[0]+x, v[1]+x, v[2]+x)

def vsubs3(v, x):
    return (v[0]-x, v[1]-x, v[2]-x)

def vadd3(v1, v2):
    return (v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2])

def vsub3(v1, v2):
    return (v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2])

def vnormalize3(v):
    m = vlength(v)
    assert(m != 0)
    return vdivs3(v, m)

def vnormalize_zero3(v):
    m = vlength(v)
    return vnormalize3(v) if m != 0 else v

#@vlogger
def vadjust_offset_length3(v, centre, min_length):
    """
    Returns a vector offset as if centre point was moved to zero; and, at least as long as min_length.
    """
    nv = vsub3(v, centre)
    m = vlength(nv)
    return nv if m > min_length else vmuls3(vnormalize_zero3(nv), min_length)

def collect_soma_nodes(pos, radius, nodes):
    """
    Creates a list of node-ids for nodes within the given soma volume.
    :param pos: centre location of soma.
    :param radius: radius of soma from pos.
    :param nodes: nodes list in skeleton data structure from amiramesh reader.
    :return: list of node-ids for nodes within soma region.
    """
    logging.info('Soma pos:%s radius:%s', pos, radius)

    soma_ids = []
    rsqr = square(radius)
    for nidx, node in nodes.iteritems():
        npos = (node.x, node.y, node.z)
        dsqr = distance_squared(pos, npos)

        if dsqr <= rsqr:
            soma_ids.append(nidx)
            logging.debug('Soma Node:%s pos:%s', nidx, npos)

    return soma_ids

def collect_node_positions(nodes):
    """
    Creates a list of node position tuples for nodes.
    :param nodes: nodes list in skeleton data structure from amiramesh reader.
    :return: list of node positions.
    """
    return [(node.x, node.y, node.z) for nidx, node in nodes.iteritems()]

def show_node_pos_stats(nodepositions):
    x = map(lambda a: a[0], nodepositions)
    y = map(lambda a: a[1], nodepositions)
    z = map(lambda a: a[2], nodepositions)

    logging.info( "X max:%s min:%s avg:%s", max(x), min(x), sum(x)/float(len(x)))
    logging.info( "Y max:%s min:%s avg:%s", max(y), min(y), sum(y)/float(len(y)))
    logging.info( "Z max:%s min:%s avg:%s", max(z), min(z), sum(z)/float(len(z)))

def debug_soma(soma, radius):
    """
    Grows fake soma nodes to outline soma visually.  Invoke prior to adding soma points.
    Assumes centre is (0,0,0)
    :param soma: BBPSDK Soma object
    :param radius: Soma radius
    """

    k_POINTS = 25

    # axis
    n = soma.grow(radius*2, 0, 0, 0.1, Section_Type.DENDRITE)
    n = soma.grow(0,radius*2, 0, 0.1, Section_Type.DENDRITE)
    n.grow(1, radius*2, 0, 0.1, Section_Type.DENDRITE)
    n = soma.grow(0,0,radius*2, 0.1, Section_Type.DENDRITE)
    n.grow(0, 1, radius*2, 0.1, Section_Type.DENDRITE)
    n.grow(1, 0, radius*2, 0.1, Section_Type.DENDRITE)

    # exterior
    for a in range(0,k_POINTS):
        ang = a * (360.0 / k_POINTS)
        i = math.sin(ang) * sradius
        j = math.cos(ang) * sradius
        n = soma.grow(i,j,0, 0.1, Section_Type.DENDRITE)
        n = soma.grow(i,0,j, 0.1, Section_Type.DENDRITE)
        n = soma.grow(0,i,j, 0.1, Section_Type.DENDRITE)


def create_node_graph(skel):
    """
    Creates a bidirectional graph dictionary of edges mapping node-id to node-ids.
    :param skel: skeleton data structure from amiramesh reader
    :return: bidirectional edge dictionary mapping node-id to set of node-ids.
    """
    edges = defaultdict(lambda: set())
    for segm in skel.segments:
        edges[segm.start].add(segm.end)
        edges[segm.end].add(segm.start)
    return edges

def create_directed_graph(somanodes, nodesgraph, allowcycles=True, somaconnected=True):
    """
    Creates a directed graph dictionary of edges mapping node id to node ids.
    :param somanodes: list of soma node-ids
    :param nodesgraph: bidirectional node-id graph of skeleton structure
    :param allowcycles: boolean set True to allow cyclic graphs, False forces acyclic graph.
    :param somaconnected: boolean set True to allow soma nodes to connect to each other, False makes them root nodes.
    :return: directed edge dictionary mapping node-id to set of node-ids.
    """
    def node_name(n, snodes, vnodes):
        return "%s%snode %s" % ('visited ' if n in vnodes else '', 'soma ' if n in snodes else '', n)

    edges = {}
    visited = []
    frontier = copy.deepcopy(somanodes)
    while frontier:
        n = frontier[0]
        neighbours = nodesgraph[n]

        logging.debug("Exploring frontier node:%s neighbours:%s", n, neighbours)

        visited.append(n)
        frontier.remove(n)
        if (n not in edges):
            edges[n] = set()

        for nn in neighbours:
            is_visited = nn in visited
            if (not is_visited):
                frontier.append(nn)
            if (somaconnected or nn not in somanodes) and (allowcycles or not is_visited):
                edges[n].add(nn)
            else:
                logging.warning("Warning - Ignoring edge from %s to %s",
                      node_name(n, somanodes, visited), node_name(nn, somanodes, visited))

    return edges

def create_node_segments_dict(segments, dgraph):
    """
    Creates a dictionary of correctly ordered segments ordered according to the dgraph.
    :param segments: list of segments from amiramesh reader.
    :param dgraph: directed node-id graph.
    :return: dictionary mapping start node-ids to the segments which grow from them.
    """
    nodesegments = defaultdict(lambda: [])
    for s in segments:
        connected = False
        if (s.start in dgraph and s.end in dgraph[s.start]):
            nodesegments[s.start].append(s)
            connected = True
        if (s.end in dgraph and s.start in dgraph[s.end]):
            r = copy.deepcopy(s)
            r.start, r.end = s.end, s.start
            r.points = [i for i in reversed(s.points)]
            nodesegments[r.start].append(r)
            connected = True

        if not connected:
            logging.warning("WARNING - unconnected segment %s->%s with %i points", s.start, s.end, len(s.points))

    return nodesegments

def validate_graph_segments(dgraph, nodesegments, somanodes = None):
    """
    Verify that the dgraph and nodesegments agree.
    :param dgraph: directed node-id graph
    :param nodesegments: dictionary mapping start node-ids to the segments which grow from them.
    :param somanodes: Optional, list of (possible) soma node-ids not to connect to

    """
    for nidx, cnodes in dgraph.iteritems():
        # nodes without children, or only somanode children, do not need segments
        assert((not cnodes or nidx in nodesegments) or
               (somanodes and all([c in somanodes for c in cnodes]))),\
                'Node %s -> [%s] from directed graph is missing from node segments dictionary.' % (nidx, cnodes)
        nsendidxs = [i.end for i in nodesegments[nidx]]
        for cidx in cnodes:
            assert(cidx in nsendidxs)

    for nidx, csegms in nodesegments.iteritems():
        assert(nidx in dgraph)
        nsstartidxs = [i.start for i in csegms]
        nsendidxs = [i.end for i in csegms]
        # if node is a leaf node, its node segments can be empty
        assert(not nsstartidxs or set(nsstartidxs) == set([nidx])), \
                'Start nodes [%s] for segments of node %s must match.' % (set([nsstartidxs]), nidx)
        for cidx in nsendidxs:
            assert(cidx in dgraph)



def grow_soma(soma, somanodes, nodesegments, nodes):
    """
    Grows the soma nodes.
    :param soma: BBPSDK Soma object
    :param somanodes: list of soma node-ids
    :param nodesegments: list of segments from amiramesh reader.
    :param nodes: dictionary mapping node positions to BBPSDK section.
    :return: directed edge dictionary mapping node-id to list of node-ids.
    """
    soma_spoints = soma.surface_points()

    # TODO: parameterize
    k_INFLATE_SOMA = False      # Default: True

    # visual debug support
    if logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
        debug_soma(soma, sradius)

    # initialize soma and nodes
    for snode_idx in somanodes:
        segments = nodesegments[snode_idx]
        for segm in segments:
            assert(segm.start == snode_idx)

            ndata = segm.points[0]
            npos = (ndata.x, ndata.y, ndata.z)
            # TODO: Fix vadjust
            snpos = vadjust_offset_length3(npos, scentre, sradius)

            if k_INFLATE_SOMA:
                soma_spoints.insert(Vector3f(snpos[0], snpos[1], snpos[2]))
                nodes[npos] = soma
            else:
                node = soma.grow(snpos[0], snpos[1], snpos[2], ndata.diameters[0], Section_Type.DENDRITE)
                node.move_point(0, Vector3f(snpos[0], snpos[1], snpos[2]))
                nodes[npos] = node
            logging.debug('Root Node: %s', segm.start)


def grow_segments(pnode_idx, dagnodes, nodesegments, nodes, visited, depth = -1):
    if (depth == 0):
        logging.warning("Warning - max depth reached for node: %i", pnode_idx)
        return

    if (pnode_idx in visited):
        return

    # TODO: parameterize
    k_SEGMENT_THRESHOLD_SQR = 0.001  # Default: 0
    k_CLIP_INSIDE_SOMA = True       # Default: True

    logging.debug('Growing:%s', str(pnode_idx))

    visited.append(pnode_idx)

    # grow sections for parent node
    segments = nodesegments[pnode_idx]
    for segm in segments:
        assert(segm.start == pnode_idx)

        if len(segm.points) < 2:
            continue

        # ndata is the parent node data (first in the segment); spt is the first section
        ndata = segm.points[0]
        npos = (ndata.x, ndata.y, ndata.z)

        logging.debug('Segment Start:%s End:%s', str(segm.end), str(segm.start))

        # start node should already exist
        assert(npos in nodes), 'Missing start node - id: %i, npos: %s' % (segm.start, npos)
        node = nodes[npos]

        # grow initial sections; first and last points belong to the start and end nodes
        # growth begins when segment exits soma
        section = None
        for pt in segm.points[1:-1]:
            pos = (pt.x, pt.y, pt.z)
            # TODO: Fix vadjust
            pos = vadjust_offset_length3(pos, scentre, 0)

            if section:
                if (distance_squared(pos, prev_pos) >= k_SEGMENT_THRESHOLD_SQR):
                    section.grow(pos[0], pos[1], pos[2], pt.diameters[0])
                    prev_pos = pos
                else:
                    logging.warning("Warning - ignoring pos: %s too close to previous: %s", pos, prev_pos)
            elif not k_CLIP_INSIDE_SOMA or vlength(pos) > sradius+pt.diameters[0]:
                section = node.grow(pos[0], pos[1], pos[2], pt.diameters[0], Section_Type.DENDRITE)
                prev_pos = pos

        # end node
        ndata = segm.points[-1]
        npos = (ndata.x, ndata.y, ndata.z)
        # TODO: Fix vadjust
        nposadj = vadjust_offset_length3(npos, scentre, max(0, sradius-ndata.diameters[0]))

        if not section and (not k_CLIP_INSIDE_SOMA or vlength(nposadj) >= sradius + ndata.diameters[0]):
            section = node

        if npos not in nodes:
            if section:
                nodes[npos] = section.grow(nposadj[0], nposadj[1], nposadj[2], ndata.diameters[0], Section_Type.DENDRITE)
            else:
                nodes[npos] = node
            logging.debug('New Node:%s', str(segm.end))

    # grow children
    for cn_idx in dagnodes[pnode_idx]:
        grow_segments(cn_idx, dagnodes, nodesegments, nodes, visited, depth - 1 if depth > 0 else -1)


def create_morphology(skel, soma_data, depth = -1):
    """
    creates morphology from the skeleton obtained
    :param skel: skeleton data structure from amiramesh reader
    :param soma_data: soma data dictionary
    :param depth: constrain the growth size to max depth size deep (or unlimited if -1)
    :return: BBPsdk morphology of the skeleton
    """

    # Collect soma nodes
    soma_center = (soma_data['centre']['x'], soma_data['centre']['y'], soma_data['centre']['z'])
    soma_radius = soma_data['radius']

    # TODO: parameterize
    k_CONNECT_SOMA_SOMA = False  # Default: False
    k_ALLOW_CYCLES = True # Default: False

    # TODO: Fix hack
    global scentre
    global sradius
    scentre = soma_center
    sradius = soma_radius
    depth = -1

    soma_node_idxs = collect_soma_nodes(soma_center, soma_radius, skel.nodes)

    npositions = collect_node_positions(skel.nodes)
    show_node_pos_stats(npositions)

    logging.info('Collected %s soma nodes out of %s total nodes',  str(len(soma_node_idxs)), str(len(skel.nodes)))

    # Create graph / data-structures of skeleton
    # NOTE: creating the directed graph also re-orders the segment directions (required to grow correctly)
    node_idx_graph = create_node_graph(skel)
    dag_nodes = create_directed_graph(soma_node_idxs, node_idx_graph, k_ALLOW_CYCLES, k_CONNECT_SOMA_SOMA)
    node_segments = create_node_segments_dict(skel.segments, dag_nodes)

    # TODO: add better tools for analysing the connectivity of unreachable nodes
    # some nodes are unreachable islands in the graph (no path from the soma); we validate and warn
    validate_graph_segments(dag_nodes, node_segments, soma_node_idxs if k_CONNECT_SOMA_SOMA else None)

    # Grow nodes
    morphology = Morphology()
    soma = morphology.soma()
    nodes = {}

    # Grow soma nodes
    grow_soma(soma, soma_node_idxs, node_segments, nodes)

    # Grow segments from inside (soma nodes) out
    visited = []
    for snidx in soma_node_idxs:
        logging.debug('Growing Soma Node:%s', str(snidx))
        grow_segments(snidx, dag_nodes, node_segments, nodes, visited, depth)

    return morphology


def create_morphology_file(morphology, filespec):
    """
    Writes the morphology object into the specified hdf5 file.
    :param morphology: BBPSDK Morphology object.
    :param filespec: Object specifying label, filepath, and output directory
    """

    # handle existing output file
    try:
        if filespec.force_overwrite:
            os.remove(filespec.skel_out_file)
    except OSError:
        pass

    morphology.label(filespec.skel_name)

    # write file to directory
    try:
        writer = Morphology_Writer()
        writer.open(filespec.skel_out_path)
        writer.write(morphology, Morphology_Repair_Stage.RAW_MORPHOLOGY)
    except OSError:
        pass


if __name__ == '__main__':

    class options:
        force_overwrite = False
        skel_path = "."
        skel_name = None
        skel_out_path = None

        skel_am_file = None
        skel_json_file = None
        skel_out_file = None

        verbosity_level = logging.ERROR
        graph_depth = -1

        @staticmethod
        def set_pathname(arg):
            options.skel_path = os.path.abspath(os.path.dirname(arg))
            if arg[-3:] == '.am':
                options.skel_name = os.path.basename(arg[:-3])
            else:
                options.skel_name = os.path.basename(arg)

        @staticmethod
        def set_filepaths():
            if not options.skel_out_path:
                options.skel_out_path = options.skel_path

            options.skel_am_file = os.path.join(options.skel_path, options.skel_name + '.am')
            options.skel_json_file = os.path.join(options.skel_path, options.skel_name + '.annotations.json')
            options.skel_out_file = os.path.join(options.skel_out_path, options.skel_name + '.h5')

        @staticmethod
        def validate():
            if not options.skel_name:
                logging.error('ERROR - Missing skeleton name.')
                sys.exit(2)
            if not os.path.exists(options.skel_am_file):
                logging.error('ERROR - Missing source file: %s', options.skel_am_file)
                sys.exit(2)
            if not os.path.exists(options.skel_json_file):
                logging.error('ERROR - Missing annotation file: %s', options.skel_json_file)
                sys.exit(3)
            if not options.force_overwrite and os.path.exists(options.skel_out_file):
                logging.error('ERROR - Existing output file (requires force overwrite): %s', options.skel_out_file)
                sys.exit(4)

    k_FORMAT = "%(message)s" # "%(asctime)-15s %(message)s"
    logging.basicConfig(format=k_FORMAT, level=options.verbosity_level)

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hfs:o:v:",["skeleton=","output_dir=","verbose="])
    except getopt.GetoptError:
        print 'skeletonize.py -h'
        sys.exit(2)
    else:
        for opt, arg in opts:
            if opt == '-h':
                print 'Skeletonize converts an Amiramesh skeleton graph, plus annotations, into a BBPSDK cell morphology.'
                print '\nUsage:'
                print ' skeletonize.py <skeleton>'
                print ' skeletonize.py -s <skeleton> [-f] [-o <output_dir>]'
                print '\t -s <filename>\t Input skeleton filename'
                print '\t -f \t\t Force overwrite of output files'
                print '\t -v <level>\t Set verbosity level: %i-%i' % (logging.NOTSET, logging.FATAL)
                print '\t -o <dirname>\t Output directory'
                print '\nExample:'
                print '\t # creates /<path>/cell.Smt.SptGraph.h5 from /<path>/cell.Smt.SptGraph'
                print '\t skeletonize.py -s cell.Smt.SptGraph'
                print '\nNotes:'
                print '\t For input source <filename>, expected input files are:'
                print '\t\t <filename>.am # Amiramesh text file of skeleton graph'
                print '\t\t <filename>.annotations.json # JSON file with {"soma": {"centre":{"x":x,"y":y,"z":z}, "radius":r}}'
                print '\t Output file(s) are:'
                print '\t\t <filename>.h5 # BBPSDK HDF5 format'
                print '\t Verbosity levels(s) are:'
                print '\t\t all=0, debug=10, info=20, warning=30, ERROR=40'
                print '\t\t ERROR is the default logging level'
                print '\t\t Debug logging levels include visual debugging artifacts added to the morphology file.'
                print '\t\t Visual debugging includes:'
                print '\t\t\t Soma star: representation of soma size and location.'
                print '\t\t\t Coordinate axis: X, Y, Z are represented as three bars with end-fingers (0=X,1=Y,2=Z).'
                print '\t Display in rtneuron-app.py using: display_morphology_file(\'/<path>/<filename>.h5\')'
                sys.exit()
            elif opt == '-f':
                options.force_overwrite = True
            elif opt in ('-v', "--verbose"):
                options.verbosity_level = int(arg)
                logging.getLogger().setLevel(options.verbosity_level)
                logging.info("Verbosity set to: %i", options.verbosity_level)
            elif opt in ("-s", "--skeleton"):
                options.set_pathname(arg)
            elif opt in ("-o", "--output_dir"):
                options.skel_out_path = arg
                if (not os.path.isdir(options.skel_out_path)):
                    print 'ERROR - Output directory must be directory:' + options.skel_out_path
                    sys.exit(4)

        if not opts:
            if len(sys.argv) != 2:
                print 'Expected skeleton source.'
                print 'skeletonize.py -h'
                sys.exit(2)
            options.set_pathname(sys.argv[1])

        options.set_filepaths()

        options.validate()

        print 'HDF5 Skeletonizer'
        print '\t Source graph: ' + options.skel_am_file
        print '\t Source annotations: ' + options.skel_json_file
        if options.force_overwrite:
            print '\nFORCING OVERWRITE of output file: ' + options.skel_out_file + '\n'

        reader = AmirameshReader()

        with open(options.skel_am_file, 'r') as f:
            skel = reader.parse(f)

        with open(options.skel_json_file, 'r') as f:
            data = json.load(f)

        morphology = create_morphology(skel, data['soma'], options.graph_depth)

        create_morphology_file(morphology, options)

        print 'Wrote out file: ' + options.skel_out_file

    finally:
        logging.shutdown()

    """
    rtneuron-app.py
    display_morphology_file('/home/holstgr/Development/Skeletonizer/oligodandrocyte/GeometrySurface.Smt.SptGraph.am.juan.h5')
    display_morphology_file('/home/holstgr/Development/Skeletonizer/oligodandrocyte/GeometrySurface.Smt.SptGraph.am.h5')
    """
