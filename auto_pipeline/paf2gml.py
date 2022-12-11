#!/usr/bin/env python

import sys
import getopt
from pathlib import Path
from igraph import Graph

def usage():
    sys.stderr.write("\nUsage: paf2gml.py [options] <in.paf>\n\n")
    sys.stderr.write("Options:\n")
    sys.stderr.write("    -o FILE    output file [stdout]\n")
    sys.stderr.flush()
    return -1

def main(argc, argv):

    if argc < 2: return usage()

    try: opts, args = getopt.gnu_getopt(argv[1:], "o:h")
    except getopt.GetoptError as err:
        sys.stderr.write("error: {}\n".format(err))
        return usage()

    output_fname = None

    for o, a in opts:
        if o == "-h": return usage()
        elif o == "-o": output_fname = a

    paf_path = Path(args[0])

    if not paf_path.is_file():
        sys.stderr.write("error: file '{}' not found\n".format(str(paf_path)))
        sys.stderr.flush()
        return -1

    if not output_fname is None: f = open(output_fname, "w")
    else: f = sys.stdout

    vertices = []
    edges = []
    vertex_set = set()

    for line in open(str(paf_path), "r"):
        tokens = line.rstrip().split()
        source_name, target_name = tokens[0], tokens[5]
        source_length, target_length = int(tokens[1]), int(tokens[6])

        if not source_name in vertex_set:
            vertex_set.add(source_name)
            vertices.append({"name" : source_name, "length" : source_length})

        if not target_name in vertex_set:
            vertex_set.add(target_name)
            vertices.append({"name" : target_name, "length" : target_length})

        source_beg, target_beg = int(tokens[2]), int(tokens[7])
        source_end, target_end = int(tokens[3]), int(tokens[8])
        strand = tokens[4]

        edges.append({"source" : source_name, "target" : target_name, "source_beg" : source_beg, "source_end" : source_end, "target_beg" : target_beg, "target_end" : target_end, "strand" : strand})

    graph = Graph.DictList(vertices, edges, directed=True, vertex_name_attr="name", edge_foreign_keys=('source', 'target'), iterative=False)

    del graph.es["source"]
    del graph.es["target"]

    graph.write_gml(f, creator="paf2gml.py")

    if not output_fname is None: f.close()

    return 0


if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
