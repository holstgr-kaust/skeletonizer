# README #

Skeletonizer is a Python tool for converting an Amiramesh skeleton graph, plus annotations, into a BBPSDK cell morphology.

## Usage ##

skeletonize.py -h
skeletonize.py <skeleton>
skeletonize.py -s <skeleton> [-f] [-o <output_dir>]

## Examples ##

Creates /<path>/cell.Smt.SptGraph.h5 from /<path>/cell.Smt.SptGraph
skeletonize.py -s cell.Smt.SptGraph

## Notes ##

For input source <filename>, expected input files are:

* <filename>.am # Amiramesh text file of skeleton graph
* <filename>.annotations.json # JSON file with {"soma": {"centre":{"x":x,"y":y,"z":z}, "radius":r}}

Output file(s) are:

* <filename>.h5 # BBPSDK HDF5 format'

Display in rtneuron-app.py using: display_morphology_file('/<path>/<filename>.h5')

