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

**Important:** The 'display_morphology_file' requires either a relative or absolute path, not just a filename.  Without a path, the morphology may appear to load, but fail to display.


### annotations.json ###

Example annotations.json file:
```
#!json
{
  "soma": {"centre":{"x":x,"y":y,"z":z}, "radius":r},
  "stack": {"AABB":{"v1":{"x":x,"y":y,"z":z}, "v2":{"x":x,"y":y,"z":z}}}
  "skeletonize": {"threshold_segment_length":t}
}
```
**Important:** Processing data from Blender and Avizo changes the coordinate system.  XYZ coordinates in Blender become -XZY in Avizo, and thus in the Skeletonizer morphology specification.  E.g., [2, 4, 8] in Blender becomes [-2, 8, 4] in Avizo, and are specified as {"x":-1,"y":8,"z":4} in the annotations file.  Manually transform positions from Blender into the Avizo coordinate system before specifying them for Skeletonizer.

*skeletonize.threshold_segment_length* specifies the edge default edge thresholding value (like **-t** command line).  It is better to specify this information as metadata to preserve the value for later reference when analyzing the resulting morphology. 

## Cross-Sectional Data ##

The approximation of cross-sectional area / perimeter from Avizo is a gross under-estimate of the minimal distance of a cross-section, which is represented in the BBPSDK as a cylinder, for cells whose cross-sections are far from circular.

`skeleton_annotate_csv.py` is a Blender script to extract accurate cross-sectional data and (for now) store it in tab-delimted `*.csv` files.  This script generates cross-sectional data from a skeletonization
representation in an Amiramesh text file (`*.am`), and a Blender project (`*.blend`) containing the
mesh for the corresponding skeletonized object.

### Notes ###

* The script depends upon the [object_cross_section](https://developer.blender.org/T34142) Blender addon which must be installed prior to running the skeleton_annotate script.
    * The `object_cross_section.py` addon script is available for installation from the `skeletonizater/addons` directory in the Skeletonizer project. 
* It is important that the Blender project file `*.blend` contains the named object; also, the object must be in the correct position (typically this is the object whose mesh was used to create the skeletonization).
* The coordinate systems differ between Blender and Avizo (the script accounts for this).
* Blender doesn't release deleted meshes, so it is better not to chunk multiple nodes (memory usage drastically grows for typical cells)
* Blender will fail if it doesn't get all the cores it expects, use `-t 1`, and don't run more copies of Blender than real cores.
    * For now, run a one node test run on the cell to get the total number of nodes (part of the output file name).
    * If there are N nodes, then use `echo $(seq 0 1 $(( N-1 )))` to iterate over the node chunks.
    * Each invocation of the script creates a tab-delimited CSV file in the same directory as the skeleton `*.am` file.
        * The file is named after the chunk of segment cross-sectional data it contains
        * Combine these files together after running all the scripts.

Below are recipes for running this script as a single invocation, and in parallel. 

Single invocation of script

```
#!bash

blender -t 1 -b $(STACK_PROJECT).blend -P skeleton_annotate_csv.py -- "$OBJECT_NAME" $(ABS_SKELETON_FILEPATH).am $START_NODE $NUMBER_OF_NODES
```

Parallel invocation of script

```
#!bash

echo $(seq 0 1 $(( $TOTAL_NODE_COUNT - 1)) ) | xargs -d " " -n 1 -P $NUM_CORES -I{} sh -c 'blender -t 1 -b $(STACK_PROJECT).blend -P skeleton_annotate_csv.py -- "$OBJECT_NAME" $(ABS_SKELETON_FILEPATH).am $(( {} )) 1'
```

Combine *.csv data files into a single file (file is named *.dat so it isn't re-read recursively; rename to *.csv file afterwards)
```
#!bash

head -n +1 "`ls *.csv | head -1`" > cross_section.csv.dat
for i in *.csv ; do tail -n +2 "$i" >> cross_section.csv.dat ; done
```


