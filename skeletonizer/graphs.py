"""
    Skeletonizer: Python Cell Morphology Analysis and Construction Toolkit

    KAUST, BESE, Neuro-Inspired Computing Project
    (c) 2014-2015. All rights reserved.
"""
"""
    Skeletonize graphs module.
"""

import os
import re
import sys
import math
import getopt
import copy
import json
import logging
import operator
from collections import defaultdict

try:
    import skeletonizer
except ImportError:
    sys.path.append(os.path.abspath(os.path.dirname(os.path.abspath(os.path.split(__file__)[0]))))

from skeletonizer.amiramesh import *
from skeletonizer.maths import *


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
        npos = node.position()
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
    return [node.position() for nidx, node in nodes.iteritems()]


def is_cut_point(pos, aabb):
    """
    Test if a node (node_idx) is on the boundary of a cut edge (aabb).
    :param pos: position vector.
    :param aabb: The AABB (Axis Aligned Bounding Box) representing the cut boundary.
    :return: True if node position represents a cut-point, or boundary node to missing data; False otherwise.
    """
    if not aabb:
        return False

    is_cut = not inside_aabb(aabb, pos)
    return is_cut


def show_node_pos_stats(nodepositions, aabb, centre):
    x = map(lambda a: a[0], nodepositions)
    y = map(lambda a: a[1], nodepositions)
    z = map(lambda a: a[2], nodepositions)

    logging.info( "X min:%s max:%s avg:%s", min(x), max(x), sum(x)/float(len(x)))
    logging.info( "Y min:%s max:%s avg:%s", min(y), max(y), sum(y)/float(len(y)))
    logging.info( "Z min:%s max:%s avg:%s", min(z), max(z), sum(z)/float(len(z)))

    clipped_nodepositions = filter(lambda v: not inside_aabb(aabb, v), nodepositions)
    logging.info( "Stack AABB clipped nodes:%s", len(clipped_nodepositions))

    for np in clipped_nodepositions:
        anp = vadjust_offset_length3(np, centre, 0)
        logging.warning( "\t clipped node pos:%s, original source pos:%s", anp, np)


def show_graph_stats(dag_nodes, node_segments):
    """
    :param dag_nodes: directed edge dictionary mapping node-id to set of node-ids.
    :param node_segments: dictionary mapping start node-ids to the segments which grow from them.
    """
    csize = 20
    ecnts = [len(ns) for _, ns in dag_nodes.iteritems()]
    ncnts = [len(ns) for _, ns in node_segments.iteritems()]
    necnts = [len(ns) for i, ns in node_segments.iteritems() if i in dag_nodes]
    len_ecnts = len(ecnts)
    len_necnts = len(necnts)
    len_ncnts = len(ncnts)

    logging.info("DAG (%s) Edge count:%s min:%s max:%s avg:%s",
                 len_ecnts, sum(ecnts),
                 min(ecnts) if len_ecnts > 0 else 'N/A',
                 max(ecnts) if len_ecnts > 0 else 'N/A',
                 sum(ecnts)/float(len_ecnts) if len_ecnts > 0 else 'N/A')
    logging.info("DAG Node Segment (%s) Edge count:%s Edge(s) min:%s max:%s avg:%s",
                 len_necnts, sum(necnts),
                 min(necnts) if len_necnts > 0 else 'N/A',
                 max(necnts) if len_necnts > 0 else 'N/A',
                 sum(necnts)/float(len_necnts) if len_necnts > 0 else 'N/A')
    logging.info("Node Segment (%s) Edge count:%s Edge(s) min:%s max:%s avg:%s",
                 len_ncnts, sum(ncnts),
                 min(ncnts) if len_ncnts > 0 else 'N/A',
                 max(ncnts) if len_ncnts > 0 else 'N/A',
                 sum(ncnts)/float(len_ncnts) if len_ncnts > 0 else 'N/A')

    posdict = defaultdict(lambda : 0)
    for _, ns in node_segments.iteritems():
        for seg in ns:
            for pt in seg.points:
                posdict[pt] += 1
    dpcnts = [cnt for i, cnt in posdict.iteritems() if cnt > 1 ]
    dpcnts.sort()
    dpcnts.reverse()

    if dpcnts:
        logging.info("Unique Segment positions:%s Duplicate positions:%s Total non-unique positions:%s Duplicates per position min:%s max:%s avg:%s",
                        len(posdict), len(dpcnts), sum(dpcnts), min(dpcnts), max(dpcnts), sum(dpcnts)/float(len(dpcnts)))
        logging.info(" Max %s position duplicate counts:%s", csize, dpcnts[:csize])
    else:
        logging.info("No duplicate positions in node and segment graphs.")

def show_grow_stats(stats, soma):
    """
    :param stats: Statistics data to log.
    :param soma: BBPSDK Soma object
    """

    node_grow_stats = stats.node_grow_stats

    csize = 30
    gcnts = [len(ns) for _, ns in node_grow_stats.iteritems()]
    gcnts.sort()
    gcnts.reverse()

    if gcnts:
        logging.info("Grown Node (%s) Edge count:%s Edge(s) min:%s max:%s avg:%s", len(gcnts), sum(gcnts), min(gcnts), max(gcnts), sum(gcnts)/float(len(gcnts)))
        if soma in node_grow_stats:
            logging.info(" Soma grown node counts:%s", len(node_grow_stats[soma]))
        else:
            logging.info("WARNING - No nodes grown from Soma")

        logging.info(" Min %i counts:%s", csize, gcnts[-csize:])
        logging.info(" Max %i counts:%s", csize, gcnts[:csize])

        bcnts = [i for i in gcnts if i > 1]

        if bcnts:
            logging.info("Grown Branching Node (%s) Edge count:%s Edge(s) min:%s max:%s avg:%s", len(bcnts), sum(bcnts), min(bcnts), max(bcnts), sum(bcnts)/float(len(bcnts)))
            logging.info(" Min %i counts:%s", csize, bcnts[-csize:])
            logging.info(" Max %i counts:%s", csize, bcnts[:csize])
        else:
            logging.warning("No Grown Branching Nodes")
    else:
        logging.warning("No Grown Nodes")

def show_warning_stats(stats):
    """
    Log warning statistics data.
    :param stats: Statistics data to log.
    """
    warnings = False
    WARN_UNCONNECTED_SEGMENTS_cnt = stats.warn_counts[stats.k_WARN_UNCONNECTED_SEGMENTS]
    WARN_IGNORED_EDGES_cnt = stats.warn_counts[stats.k_WARN_IGNORED_EDGES]
    WARN_MAX_GROW_DEPTH_REACHED_cnt = stats.warn_counts[stats.k_WARN_MAX_GROW_DEPTH_REACHED]
    WARN_CUT_NOTES_FOUND_cnt = stats.warn_counts[stats.k_WARN_CUT_NODES_FOUND]
    INFO_IGNORED_POSITIONS_cnt = stats.warn_counts[stats.k_INFO_IGNORED_POSITIONS]

    if WARN_UNCONNECTED_SEGMENTS_cnt > 0:
        warnings = True
        logging.warning("WARNING - %s Unconnected Segments (edge-islands not reachable from soma)", WARN_UNCONNECTED_SEGMENTS_cnt)

    if WARN_IGNORED_EDGES_cnt > 0:
        warnings = True
        logging.warning("WARNING - %s Ignored Edges (possible segments returning into soma or cycles in original graph)", WARN_IGNORED_EDGES_cnt)

    if WARN_MAX_GROW_DEPTH_REACHED_cnt > 0:
        warnings = True
        logging.warning("WARNING - Reached maximum grow depth %s times", WARN_MAX_GROW_DEPTH_REACHED_cnt)

    if WARN_CUT_NOTES_FOUND_cnt > 0:
        warnings = True
        logging.warning("WARNING - Found %s cut nodes", WARN_CUT_NOTES_FOUND_cnt)

    if INFO_IGNORED_POSITIONS_cnt > 0:
        warnings = True
        logging.info("INFO - %s Ignored Segments Positions (thresholded)", INFO_IGNORED_POSITIONS_cnt)

    if warnings and logging.getLogger().getEffectiveLevel() > logging.DEBUG:
        logging.warning("NOTE: To view warning and info details, enable DEBUG verbosity: -v 10")


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

def create_directed_graph(somanodes, nodesgraph, options, stats):
    """
    Creates a directed graph dictionary of edges mapping node id to node ids.
    :param somanodes: list of soma node-ids
    :param nodesgraph: bidirectional node-id graph of skeleton structure
    :param options: struct of graph options.
    :param stats: statistic collection object
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
            if (options.k_CONNECT_SOMA_SOMA or nn not in somanodes) and (options.k_ALLOW_CYCLES or not is_visited):
                edges[n].add(nn)
            else:
                stats.warn_counts[stats.k_WARN_IGNORED_EDGES] += 1
                logging.debug("WARNING - Ignoring edge from %s to %s",
                      node_name(n, somanodes, visited), node_name(nn, somanodes, visited))

    return edges

def create_node_segments_dict(segments, dgraph, stats):
    """
    Creates a dictionary of correctly ordered segments ordered according to the dgraph.
    :param segments: list of segments from amiramesh reader.
    :param dgraph: directed node-id graph.
    :param stats: statistic collection object
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
            stats.warn_counts[stats.k_WARN_UNCONNECTED_SEGMENTS] += 1
            logging.debug("WARNING - unconnected segment %s->%s with %i points", s.start, s.end, len(s.points))

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

