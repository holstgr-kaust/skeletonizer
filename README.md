# README #

Skeletonizer is a Python tool for converting an Amiramesh skeleton graph, plus annotations, into a BBPSDK cell morphology.

## Usage ##


```
#!python

skeletonize.py -h
skeletonize.py <skeleton>
skeletonize.py -s <skeleton> [-f] [-o <output_dir>] [-v <level>] [-t <threshold>]

```

## Examples ##

Creates */<path>/cell.Smt.SptGraph.h5* from */<path>/cell.Smt.SptGraph*

```
#!python


skeletonize.py -s cell.Smt.SptGraph
```


## Notes ##

For input source <filename>, expected input files are:

* <filename>.am # Amiramesh text file of skeleton graph
* <filename>.annotations.json # JSON file with {"soma": {"centre":{"x":x,"y":y,"z":z}, "radius":r}}

Output file(s) are:

* <filename>.h5 # BBPSDK HDF5 format'

Verbosity levels(s) are: all=0, debug=10, INFO=20, warning=30, error=40

* INFO is the default logging level
* Debug logging levels include visual debugging artifacts added to the morphology file:
    * Soma star: representation of soma size and location.
    * Coordinate axis: X, Y, Z are represented as three bars with end-fingers (0=X,1=Y,2=Z).
* All logging level includes additional visual debugging artifacts:
	* Soma dendrites: visual representation of original source soma skeleton.

Threshold currently specifies the minimum segment section length.

Display in rtneuron-app.py using: display_morphology_file('/<path>/<filename>.h5')
