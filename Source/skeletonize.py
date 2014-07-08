#!/usr/bin/env python

"""
This program creates a BBPSDK HDF5 morphology file from a skeletonization representation in an Amiramesh text file.
"""
import os
import re
import sys
import getopt
import copy
import json
from collections import defaultdict

from bbp_import_module import *
from amiramesh import *


def square(x):
    return x * x

def distance_squared(p1, p2):
    return sum(map(lambda x, y: square(x - y), p1, p2))


def collect_soma_nodes(pos, radius, nodes):
    """
    Creates a list of node-ids for nodes within the given soma volume.
    :param pos: centre location of soma.
    :param radius: radius of soma from pos.
    :param nodes: nodes list in skeleton data structure from amiramesh reader.
    :return: list of node-ids for nodes within soma region.
    """
    print('Soma pos:' + str(pos) + ' radius:' + str(radius))
    soma_ids = []
    rsqr = square(radius)
    for nidx, node in nodes.iteritems():
        npos = (node.x, node.y, node.z)
        dsqr = distance_squared(pos, npos)

        if dsqr <= rsqr:
            soma_ids.append(nidx)
            print('Soma Node:' + str(nidx) + ' pos:' + str(npos))

    return soma_ids

def collect_node_positions(nodes):
    """
    Creates a list of node position tuples for nodes.
    :param nodes: nodes list in skeleton data structure from amiramesh reader.
    :return: list of node positions.
    """
    return [(node.x, node.y, node.z) for nidx, node in nodes.iteritems()]

def print_node_pos_stats(nodepositions):
    x = map(lambda a: a[0], nodepositions)
    y = map(lambda a: a[1], nodepositions)
    z = map(lambda a: a[2], nodepositions)

    print "X max:", max(x), "min:", min(x), "avg:", sum(x)/float(len(x))
    print "Y max:", max(y), "min:", min(y), "avg:", sum(y)/float(len(y))
    print "Z max:", max(z), "min:", min(z), "avg:", sum(z)/float(len(z))


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
    edges = {}
    visited = []
    frontier = copy.deepcopy(somanodes)
    while frontier:
        n = frontier[0]
        neighbours = nodesgraph[n]

        print("Exploring frontier node: " + str(n) + " neighbours:" + str(neighbours))

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
                print("Warning - Ignoring edge to node: " + str(nn))

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
            print("WARNING - unconnected segment " + str(s.start) + "->" + str(s.end) + " with " + str(len(s.points)) + " points")
            # raise Exception('segment connecting nodes does not exist in directed connection graph')

    return nodesegments

def validate_graph_segments(dgraph, nodesegments):
    """
    Verify that the dgraph and nodesegments agree
    """
    for nidx, cnodes in dgraph.iteritems():
        assert(nidx in nodesegments)
        nsendidxs = [i.end for i in nodesegments[nidx]]
        for cidx in cnodes:
            assert(cidx in nsendidxs)

    for nidx, csegms in nodesegments.iteritems():
        assert(nidx in dgraph)
        nsstartidxs = [i.start for i in csegms]
        nsendidxs = [i.end for i in csegms]
        assert(set(nsstartidxs) == set([nidx]))
        for cidx in nsendidxs:
            assert(cidx in dgraph)


def grow_soma(soma, somanodes, segments, nodes):
    """
    Grows the soma nodes.
    :param soma: BBPSDK Soma object
    :param somanodes: list of soma node-ids
    :param segments: list of segments from amiramesh reader.
    :param nodes: dictionary mapping node positions to BBPSDK section.
    :return: directed edge dictionary mapping node-id to list of node-ids.
    """
    for segm in segments:
        if segm.start not in somanodes:
            continue

        ndata = segm.points[0]
        npos = (ndata.x, ndata.y, ndata.z)
        # TODO: specify these nodes are part of the soma
        #node = soma.grow(npos[0], npos[1], npos[2], ndata.diameters[0], Section_Type.SOMA)
        node = soma.grow(npos[0], npos[1], npos[2], ndata.diameters[0], Section_Type.DENDRITE)
        # TODO: these nodes are within the soma; do we need to move them?
        node.move_point(0, Vector3f(npos[0], npos[1], npos[2]))
        nodes[npos] = node
        print('Root Node:' + str(segm.start))


def grow_segments(pnode_idx, dagnodes, nodesegments, nodes, visited, depth = -1):
    if (depth == 0):
        print("Warning - max depth reach for node: " + str(pnode_idx))
        return

    if (pnode_idx in visited):
        return

    #print('Growing:'+str(pnode_idx))

    visited.append(pnode_idx)

    # grow sections for parent node
    segments = nodesegments[pnode_idx]
    for segm in segments:
        assert(segm.start == pnode_idx)

        if len(segm.points) < 2:
            continue

        # ndata is the parent node data (first in the segment); spt is the first section
        ndata, spt = segm.points[0:2]
        npos = (ndata.x, ndata.y, ndata.z)

        #print('Segment Start:'+str(segm.start)+' End:'+str(segm.end))

        node = nodes[npos]

        section = node.grow(spt.x, spt.y, spt.z, spt.diameters[0], Section_Type.DENDRITE)

        # TODO: Sections should end at a node (unless there are loops)
        for pt in segm.points[1:]:
            pos = (pt.x, pt.y, pt.z)
            if (pos in nodes):
                section.grow(pt.x, pt.y, pt.z, pt.diameters[0])
                # TODO: if (pos in nodes) section.grow(nodes[pos]) else
                #section.grow(nodes[pos])
            else:
                section.grow(pt.x, pt.y, pt.z, pt.diameters[0])

        ndata = segm.points[-1:][0]
        npos = (ndata.x, ndata.y, ndata.z)
        if npos not in nodes:
            nodes[npos] = section.grow(npos[0], npos[1], npos[2], ndata.diameters[0], Section_Type.DENDRITE)
            print('New Node:'+str(segm.end))

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

    soma_node_idxs = collect_soma_nodes(soma_center, soma_radius, skel.nodes)

    npositions = collect_node_positions(skel.nodes)
    print_node_pos_stats(npositions)

    print('Collected ' + str(len(soma_node_idxs)) + ' soma nodes out of ' + str(len(skel.nodes)) + ' total nodes')

    # Create graph / data-structures of skeleton
    # NOTE: creating the directed graph also re-orders the segment directions (required to grow correctly)
    node_idx_graph = create_node_graph(skel)
    dag_nodes = create_directed_graph(soma_node_idxs, node_idx_graph)
    node_segments = create_node_segments_dict(skel.segments, dag_nodes)

    # TODO: add better tools for analysing the connectivity of unreachable nodes
    # some nodes are unreachable islands in the graph (no path from the soma); we validate and warn
    validate_graph_segments(dag_nodes, node_segments)

    # Grow nodes
    morphology = Morphology()
    soma = morphology.soma()
    nodes = {}

    # Grow soma nodes
    grow_soma(soma, soma_node_idxs, skel.segments, nodes)

    # Grow segments from inside (soma nodes) out
    visited = []
    for snidx in soma_node_idxs:
        #print('Growing Soma Node:'+str(snidx))
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
                print 'ERROR - Missing skeleton name.'
                sys.exit(2)
            if not os.path.exists(options.skel_am_file):
                print 'ERROR - Missing source file: ' + options.skel_am_file
                sys.exit(2)
            if not os.path.exists(options.skel_json_file):
                print 'ERROR - Missing annotation file: ' + options.skel_json_file
                sys.exit(3)
            if not options.force_overwrite and os.path.exists(options.skel_out_file):
                print 'ERROR - Existing output file (requires force overwrite): ' + options.skel_out_file
                sys.exit(4)


    try:
        opts, args = getopt.getopt(sys.argv[1:],"hfs:o:",["skeleton=","output_dir="])
    except getopt.GetoptError:
        print 'skeletonize.py -h'
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print 'Skeletonize converts an Amiramesh skeleton graph, plus annotations, into a BBPSDK cell morphology.'
            print '\n Usage:'
            print 'skeletonize.py <skeleton>'
            print 'skeletonize.py -s <skeleton> [-f] [-o <output_dir>]'
            print '\t -s <filename>\t Input skeleton filename'
            print '\t -f \t\t\t Force overwrite of output files'
            print '\t -o <dirname>\t Output directory'
            print '\n Example:'
            print '\t # creates /<path>/cell.Smt.SptGraph.h5 from /<path>/cell.Smt.SptGraph'
            print '\t skeletonize.py -s cell.Smt.SptGraph'
            print '\n Notes:'
            print '\t For input source <filename>, expected input files are:'
            print '\t\t <filename>.am # Amiramesh text file of skeleton graph'
            print '\t\t <filename>.annotations.json # JSON file with {"soma": {"centre":{"x":x,"y":y,"z":z}, "radius":r}}'
            print '\t Output file(s) are:'
            print '\t\t <filename>.h5 # BBPSDK HDF5 format'
            print '\t Display in rtneuron-app.py using: display_morphology_file(\'/<path>/<filename>.h5\')'
            sys.exit()
        elif opt == '-f':
            options.force_overwrite = True
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

    morphology = create_morphology(skel, data['soma'])

    create_morphology_file(morphology, options)

    print 'Wrote out file: ' + options.skel_out_file

    """
    rtneuron-app.py
    display_morphology_file('/home/holstgr/Development/Skeletonizer/oligodandrocyte/GeometrySurface.Smt.SptGraph.am.juan.h5')
    display_morphology_file('/home/holstgr/Development/Skeletonizer/oligodandrocyte/GeometrySurface.Smt.SptGraph.am.h5')
    """
