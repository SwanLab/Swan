import struct
import sys

import numpy as np

HEADER_BYTES = 32
NODE_RECORD_BYTES = 56
NODE_COORD_OFFSET = 32
ELEMENT_HEADER_BYTES = 28
ELEMENT_RECORD_BYTES = 48
ELEMENT_NODES = slice(8, 12)
METERS_TO_MM = 1000


def read_nodes(data):
    n_nodes = struct.unpack_from('<i', data, 28)[0]
    ids, coords = [], []
    for k in range(n_nodes):
        offset = HEADER_BYTES + k*NODE_RECORD_BYTES
        ids.append(struct.unpack_from('<i', data, offset)[0])
        coords.append(struct.unpack_from('<2d', data, offset + NODE_COORD_OFFSET))
    return np.array(ids), METERS_TO_MM*np.array(coords), HEADER_BYTES + n_nodes*NODE_RECORD_BYTES


def read_quads(data, offset, node_ids):
    n_elements = struct.unpack_from('<i', data, offset + 20)[0]
    start = offset + ELEMENT_HEADER_BYTES
    quads = np.array([struct.unpack_from('<12i', data, start + k*ELEMENT_RECORD_BYTES)[ELEMENT_NODES]
                      for k in range(n_elements)])
    index = {node_id: k for k, node_id in enumerate(node_ids)}
    return np.vectorize(index.get)(quads)


def find_boundary_edges(quads):
    edges = np.concatenate([quads[:, [k, (k + 1) % 4]] for k in range(4)])
    keys, counts = np.unique(np.sort(edges, axis=1), axis=0, return_counts=True)
    single = {tuple(key) for key in keys[counts == 1]}
    return [tuple(e) for e in edges if tuple(sorted(e)) in single]


def chain_loops(edges):
    next_node = dict(edges)
    loops = []
    while next_node:
        start, node = next(iter(next_node.items()))
        loop = [start]
        del next_node[start]
        while node != start:
            loop.append(node)
            node = next_node.pop(node)
        loops.append(loop)
    return loops


def remove_collinear(points, tol=1e-9):
    previous = np.roll(points, 1, axis=0)
    following = np.roll(points, -1, axis=0)
    cross = (points - previous)[:, 0]*(following - points)[:, 1] - (points - previous)[:, 1]*(following - points)[:, 0]
    return points[np.abs(cross) > tol]


def write_loops(file_name, loops, coords):
    with open(file_name, 'w') as f:
        for loop in loops:
            for x, y in remove_collinear(coords[loop]):
                f.write('%.9f,%.9f\n' % (x, y))
            f.write('NaN,NaN\n')


def main(msh_file, out_file):
    data = open(msh_file, 'rb').read()
    node_ids, coords, offset = read_nodes(data)
    quads = read_quads(data, offset, node_ids)
    loops = chain_loops(find_boundary_edges(quads))
    write_loops(out_file, loops, coords)
    print('nodes', len(node_ids), 'quads', len(quads), 'loops', [len(loop) for loop in loops])


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
