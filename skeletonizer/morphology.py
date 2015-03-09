"""
    Skeletonizer: Python Cell Morphology Analysis and Construction Toolkit

    KAUST, BESE, Neuro-Inspired Computing Project
    (c) 2014-2015. All rights reserved.
"""
"""
    Skeletonize morphology module.
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

from skeletonizer.bbp_import_module import *
from skeletonizer.maths import *
from skeletonizer.graphs import *

class MorphologyCreateOptions:
    force_overwrite = False
    skel_path = "."
    skel_name = None
    skel_out_path = None

    skel_am_file = None
    skel_json_file = None
    skel_out_file = None

    verbosity_level = logging.INFO
    force_segment_threshold = False
    threshold_segment_length = 0
    scaling_factor = 1
    allow_cycles = False
    graph_depth = -1

    stack_AABB = None


    def set_pathname(self, arg):
        self.skel_path = os.path.abspath(os.path.dirname(arg))
        if arg[-3:] == '.am':
            self.skel_name = os.path.basename(arg[:-3])
        else:
            self.skel_name = os.path.basename(arg)

    def set_filepaths(self):
        if not self.skel_out_path:
            self.skel_out_path = self.skel_path

        self.skel_am_file = os.path.join(self.skel_path, self.skel_name + '.am')
        self.skel_json_file = os.path.join(self.skel_path, self.skel_name + '.annotations.json')
        self.skel_out_file = os.path.join(self.skel_out_path, self.skel_name + '.h5')

    def set_annotation_data(self, data):
        if 'skeletonize' in data:
            skeletonize_config = data['skeletonize']
            assert (type(skeletonize_config) == dict), \
                    "Expected skeletonize section dictionary object"
            if 'threshold_segment_length' in skeletonize_config and not self.force_segment_threshold:
                self.threshold_segment_length = float(skeletonize_config['threshold_segment_length'])
                logging.info("Segment length threshold set to: %f", self.threshold_segment_length)

        if 'stack' in data:
            stack_metadata = data['stack']
            assert (type(stack_metadata) == dict), \
                    "Expected stack section dictionary object"
            if 'AABB' in stack_metadata:
                # TODO: fix constant
                adjust_amt = -1.0
                v1_x = stack_metadata['AABB']['v1']['x']
                v1_y = stack_metadata['AABB']['v1']['y']
                v1_z = stack_metadata['AABB']['v1']['z']
                v2_x = stack_metadata['AABB']['v2']['x']
                v2_y = stack_metadata['AABB']['v2']['y']
                v2_z = stack_metadata['AABB']['v2']['z']

                self.stack_AABB = adjust_aabb(
                    v3_to_aabb((v1_x, v1_y, v1_z), (v2_x, v2_y, v2_z)), adjust_amt)

                logging.info("Found stack AABB: %f",
                             self.threshold_segment_length)
    #TODO: throw exception instead of sys.exit (client should sys.exit)
    def validate(self):
        if not self.skel_name:
            logging.error('ERROR - Missing skeleton name.')
            sys.exit(2)
        if not os.path.exists(self.skel_am_file):
            logging.error('ERROR - Missing source file: %s', self.skel_am_file)
            sys.exit(2)
        if not os.path.exists(self.skel_json_file):
            logging.error('ERROR - Missing annotation file: %s', self.skel_json_file)
            sys.exit(3)
        if not self.force_overwrite and os.path.exists(self.skel_out_file):
            logging.error('ERROR - Existing output file (requires force overwrite): %s', self.skel_out_file)
            sys.exit(4)


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
        i = math.sin(ang) * radius
        j = math.cos(ang) * radius
        n = soma.grow(i,j,0, 0.1, Section_Type.DENDRITE)
        n = soma.grow(i,0,j, 0.1, Section_Type.DENDRITE)
        n = soma.grow(0,i,j, 0.1, Section_Type.DENDRITE)

def debug_scale_cut_point_diameter(scaled_diameter, scale):
    """
    Returns a new scaled diameter for visual debugging of cut-point nodes.
    :param scaled_diameter: The current, scaled node diameter
    :param scale: The scaling factor
    :return: New diameter for cut node
    """
    return max(2 * scale, scaled_diameter * 5)


def grow_soma(soma, somanodes, nodesegments, nodes, offsets, options, stats):
    """
    Grows the soma nodes.
    :param soma: BBPSDK Soma object.
    :param somanodes: list of soma node-ids.
    :param nodesegments: dictionary mapping start node-ids to the segments which grow from them.
    :param nodes: dictionary mapping node positions to BBPSDK section.
    :param offsets: tuple of (soma_centre, soma_radius).
    :param options: struct of growth options.
    :param stats: statistic collection object
    """

    # NOTE: we offset the original graph to centre the soma at origin in BBPSDK morphology, but preserve
    # the original positions to make it easier to report original graph positions to user
    scentre, sradius = offsets
    scale = options.k_SCALING_FACTOR
    soma_spoints = soma.surface_points()

    # visual debug support
    if logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
        debug_soma(soma, sradius * scale)

    # initialize soma and nodes
    for snode_idx in somanodes:
        segments = nodesegments[snode_idx]
        for segm in segments:
            assert(segm.start == snode_idx)

            ndata = segm.points[0]
            npos = ndata.position()
            snpos = vadjust_offset_length3(npos, scentre, sradius)

            if options.k_INFLATE_SOMA:
                spos = vmuls3(snpos, scale)
                sdiameter = ndata.diameter * scale
                soma_spoints.insert(Vector3f(spos[0], spos[1], spos[2]))
                nodes[npos] = soma
            else:
                if logging.getLogger().getEffectiveLevel() < logging.DEBUG:
                    snpos = vadjust_offset_length3(npos, scentre, 0)
                spos = vmuls3(snpos, scale)
                sdiameter = ndata.diameter * scale
                node = soma.grow(spos[0], spos[1], spos[2], sdiameter, Section_Type.DENDRITE)
                stats.node_grow_stats[soma].append(snpos)
                node.move_point(0, Vector3f(spos[0], spos[1], spos[2]))
                nodes[npos] = node

            logging.debug('Root Node: %s', segm.start)

    # debug support
    logging.info("Soma created: radius (mean, max): (%s, %s)",
                    soma.mean_radius(), soma.max_radius())


def grow_segments(pnode_idx, dagnodes, nodesegments, nodes, visited,
                  morphology, offsets, options, stats, depth = -1):
    """
    Grows the node to node segments.
    :param pnode_idx: node-id of parent node.
    :param dagnodes: directed edge dictionary mapping node-id to list of node-ids.
    :param nodesegments: dictionary mapping start node-ids to the segments which grow from them.
    :param nodes: dictionary mapping node positions to BBPSDK section.
    :param visited: node-ids of already visited nodes.
    :param morphology: BBPSDK Morphology object.
    :param offsets: tuple of (soma_centre, soma_radius).
    :param options: struct of growth options.
    :param stats: statistic collection object
    :param depth: debugging: controls growth size; if non-negative specifies max node count; -1 if unlimited
    """
    if (depth == 0):
        stats.warn_counts[stats.k_WARN_MAX_GROW_DEPTH_REACHED] += 1
        logging.debug("WARNING - max depth reached for node: %i", pnode_idx)
        return

    if (pnode_idx in visited):
        return

    # NOTE: we offset the original graph to centre the soma at origin in BBPSDK morphology, but preserve
    # the original positions to make it easier to report original graph positions to user
    scentre, sradius = offsets
    scale = options.k_SCALING_FACTOR

    logging.debug('Growing:%s', str(pnode_idx))

    visited.append(pnode_idx)

    is_parent_cut = False

    # grow sections for parent node
    segments = nodesegments[pnode_idx]
    for segm in segments:
        assert(segm.start == pnode_idx)

        if len(segm.points) < 2:
            continue

        # ndata is the parent node data (first in the segment); spt is the first section
        ndata = segm.points[0]
        npos = ndata.position()

        logging.debug('Segment Start:%s End:%s', str(segm.start), str(segm.end))

        # start node should already exist
        assert(npos in nodes), 'Missing start node - id: %i, npos: %s' % (segm.start, npos)
        node = nodes[npos]

        is_parent_cut = is_parent_cut or is_cut_point(npos, options.k_CUTPOINT_AABB)
        if is_parent_cut:
            logging.debug('Cut node reached at node:%s position:%s', str(segm.start), npos)
            break

        is_cut = is_parent_cut

        # grow initial sections; first and last points belong to the start and end nodes
        # growth begins when segment exits soma
        section = None
        for pt in segm.points[1:-1]:
            pos_orig = pt.position()
            pos = vadjust_offset_length3(pos_orig, scentre, 0)

            spos = vmuls3(pos, scale)
            sdiameter = pt.diameter * scale

            is_cut = is_cut or is_cut_point(pos_orig, options.k_CUTPOINT_AABB)

            # visual debug support
            if is_cut and logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
                sdiameter = debug_scale_cut_point_diameter(sdiameter, scale)

            if section:
                if (distance_squared(pos, prev_pos) >= options.k_SEGMENT_THRESHOLD_SQR):
                    section.grow(spos[0], spos[1], spos[2], sdiameter)

                    if is_cut:
                        morphology.mark_cut_point(section)

                    prev_pos = pos
                else:
                    stats.warn_counts[stats.k_INFO_IGNORED_POSITIONS] += 1
                    logging.debug("INFO - ignoring pos: %s too close to previous: %s", pos, prev_pos)
            elif not options.k_CLIP_INSIDE_SOMA or vlength(pos) > sradius+pt.diameter:
                section = node.grow(spos[0], spos[1], spos[2], sdiameter, Section_Type.DENDRITE)
                stats.node_grow_stats[node].append(pos)
                prev_pos = pos

            if is_cut:
                stats.warn_counts[stats.k_WARN_CUT_NODES_FOUND] += 1
                logging.debug('Cut node reached at node segment position:%s', pos_orig)
                break

        # end node
        ndata = segm.points[-1]
        npos = ndata.position()
        nposadj = vadjust_offset_length3(npos, scentre, max(0, sradius-ndata.diameter))

        is_cut = is_cut or is_cut_point(npos, options.k_CUTPOINT_AABB)

        if is_cut:
            stats.warn_counts[stats.k_WARN_CUT_NODES_FOUND] += 1
            logging.debug('Ending cut node reached at node position:%s', npos)

        if not section and (not options.k_CLIP_INSIDE_SOMA or vlength(nposadj) >= sradius + ndata.diameter):
            section = node

        if npos not in nodes:
            if section:
                spos = vmuls3(nposadj, scale)
                sdiameter = ndata.diameter * scale

                # visual debug support
                if is_cut and logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
                    sdiameter = debug_scale_cut_point_diameter(sdiameter, scale)

                nodes[npos] = section.grow(spos[0], spos[1], spos[2], sdiameter, Section_Type.DENDRITE)

                if is_cut:
                    morphology.mark_cut_point(nodes[npos])

                stats.node_grow_stats[section].append(nposadj)
                logging.debug('New Node:%s', str(segm.end))
            else:
                nodes[npos] = node
                logging.debug('Reusing Start Node:%s as End Node:%s', str(segm.start), str(segm.end))

    if is_parent_cut:
        return

    # grow children
    for cn_idx in dagnodes[pnode_idx]:
        grow_segments(cn_idx, dagnodes, nodesegments, nodes, visited, morphology,
                      offsets, options, stats, depth - 1 if depth > 0 else -1)


def create_morphology(skel, soma_data, options):
    """
    creates morphology from the skeleton obtained
    :param skel: skeleton data structure from amiramesh reader
    :param soma_data: soma data dictionary
    :param options: struct of create morphology options
    :return: BBPsdk morphology of the skeleton
    """
    class morph_options:
        # boolean set True to allow cyclic graphs, False forces acyclic graph.
        k_ALLOW_CYCLES = options.allow_cycles                                           # Default: False
        # boolean set True to allow soma nodes to connect to each other, False makes them root nodes.
        k_CONNECT_SOMA_SOMA = options.verbosity_level <= logging.NOTSET                 # Default: False

        # float specifies minimum length between segment arcs (inter-node section edges)
        k_SEGMENT_THRESHOLD_SQR = options.threshold_segment_length                      # Default: 0
        # boolean set True to start segments after they leave the soma
        k_CLIP_INSIDE_SOMA = options.verbosity_level > logging.NOTSET                   # Default: True

        # boolean set True to create normal BBPSDK soma node; False creates zero sized node for debugging
        k_INFLATE_SOMA = options.verbosity_level > logging.NOTSET                       # Default: True

        # float specifies morphology scaling factor
        k_SCALING_FACTOR = options.scaling_factor                                       # Default: 1

        k_CUTPOINT_AABB = options.stack_AABB                                            # Default: None


    class morph_statistics:
        k_WARN_UNCONNECTED_SEGMENTS = 1
        k_WARN_IGNORED_EDGES = 2
        k_WARN_MAX_GROW_DEPTH_REACHED = 3
        k_WARN_CUT_NODES_FOUND = 4
        k_INFO_IGNORED_POSITIONS = 100

        # dictionary mapping warning and info categories above to occurrence counts
        warn_counts = defaultdict(lambda: 0)

        # dictionary mapping BBPSDK nodes to the positions grown from them.
        node_grow_stats = defaultdict(lambda: [])


    depth = options.graph_depth

    # Collect soma nodes
    soma_centre = (soma_data['centre']['x'], soma_data['centre']['y'], soma_data['centre']['z'])
    soma_radius = soma_data['radius']

    soma_node_idxs = collect_soma_nodes(soma_centre, soma_radius, skel.nodes)

    npositions = collect_node_positions(skel.nodes)

    show_node_pos_stats(npositions, options.stack_AABB, soma_centre)
    logging.info('Collected %s soma nodes out of %s total nodes',  str(len(soma_node_idxs)), str(len(skel.nodes)))

    # Create graph / data-structures of skeleton
    # NOTE: creating the directed graph also re-orders the segment directions (required to grow correctly)
    node_idx_graph = create_node_graph(skel)
    dag_nodes = create_directed_graph(soma_node_idxs, node_idx_graph,
                                      morph_options, morph_statistics)
    node_segments = create_node_segments_dict(skel.segments, dag_nodes,
                                              morph_statistics)

    show_graph_stats(dag_nodes, node_segments)

    # TODO: add better tools for analysing the connectivity of unreachable nodes
    # some nodes are unreachable islands in the graph (no path from the soma); we validate and warn
    validate_graph_segments(dag_nodes, node_segments,
                            soma_node_idxs if morph_options.k_CONNECT_SOMA_SOMA else None)

    # Grow nodes
    morphology = Morphology()
    soma = morphology.soma()
    nodes = {}

    # Grow soma nodes
    grow_soma(soma, soma_node_idxs,
              node_segments, nodes,
              (soma_centre, soma_radius),
              morph_options, morph_statistics)

    # Grow segments from inside (soma nodes) out
    visited = []
    for snidx in soma_node_idxs:
        logging.debug('Growing Soma Node:%s', str(snidx))
        grow_segments(snidx, dag_nodes, node_segments, nodes, visited,
                      morphology, (soma_centre, soma_radius),
                      morph_options, morph_statistics, depth)

    show_grow_stats(morph_statistics, soma)
    show_warning_stats(morph_statistics)

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

