#!/usr/bin/env python

import numpy as np
import sys
from xml.etree import ElementTree
from xml.dom import minidom

class unitcell:
    """
    Extract unitcell structure from AkaiKKR output

    Attributes
    ----------
    dimension : int
        dimension of space
    bravais : string
        type of Bravais lattice
    lattice_constant : array of floats
        lattice constans of Bravais lattice [a,b,c]
    lattice_angle : array of floats
        lattice angles of Bravais lattice [alpha,beta,gamma]
    cell_bases : dimension x dimension matrix of floats
        (column) basis vectors of unit cell
    site_labels : list of strings
        list of site labels
    sites : list of [array of floats, int]
        relative site coordinates and type of each site
    coupling_constants = list of floats
        list of coupling constants, J
    bonds : list of [int, int, array of int, int]
        source, target, cell offset of target, and type of each bond
    """
    def __init__(self, path, cutoff):
        """
        Paramters
        ---------
        path :
           path to the input file
        """
        self.parse(path, cutoff)
        self.volume = abs(np.linalg.det(self.cell_bases))

    def parse(self, path, cutoff):
        """
        Read effective lattice model information from AkaiKKR output file

        Parameters
        ----------
        path : path to the input file
        """
        self.dimension = 3
        self.bravais = ''
        self.lattice_constant = []
        self.lattice_angle = []
        self.cell_bases = np.zeros([3,3], order='F')
        self.site_labels = []
        self.sites = []
        self.coupling_constants = []
        self.bonds = []
        meVinK = 1.16045221e+1 # meV in unit of Kelvin
        section = ''
        site_type_map = {}
        with open(path) as f:
            for line in f:
                words = line.split()
                if (words == ['lattice', 'constant']):
                    section = 'lattice'
                if (section == 'lattice'):
                    if (len(words) > 3 and words[3] == 'c/a='):
                        self.bravais = words[0].split('=')[1]
                        a = float(words[2])
                        self.lattice_constant = np.array([a, a * float(words[6]), a * float(words[4])])
                    if (len(words) > 0 and words[0] == 'alpha='):
                        self.lattice_angle = np.array([float(words[1]), float(words[3]), float(words[5])])

                if (words[0:3] == ['primitive', 'translation', 'vectors']):
                    section = 'primitive'
                if (section == 'primitive'):
                    if (len(words) > 0 and words[0] == 'a=('):
                        self.cell_bases[0][0] = float(words[1])
                        self.cell_bases[1][0] = float(words[2])
                        self.cell_bases[2][0] = float(words[3].split(')')[0])
                    if (len(words) > 0 and words[0] == 'b=('):
                        self.cell_bases[0][1] = float(words[1])
                        self.cell_bases[1][1] = float(words[2])
                        self.cell_bases[2][1] = float(words[3].split(')')[0])
                    if (len(words) > 0 and words[0] == 'c=('):
                        self.cell_bases[0][2] = float(words[1])
                        self.cell_bases[1][2] = float(words[2])
                        self.cell_bases[2][2] = float(words[3].split(')')[0])

                if (words == ['type', 'of', 'site']):
                    section = 'type'
                if (section == 'type'):
                    if (len(words) > 1 and words[1] == 'rmt='):
                        t = words[0].split('=')[1]
                        i = len(self.site_labels)
                        self.site_labels.append(t)
                        site_type_map[t] = i

                if (words == ['atoms', 'in', 'the', 'unit', 'cell']):
                    section = 'atoms'
                if (section == 'atoms'):
                    if (len(words) > 0 and words[0] == 'position='):
                        self.sites.append([np.array([float(words[1]), float(words[2]), float(words[3])]),
                                           site_type_map[words[4].split('=')[1]]])

                if (words == ['J_ij']):
                    section = 'J_ij'
                if (section == 'J_ij'):
                    if (len(words) == 12):
                        self.coupling_constants.append(2 * float(words[10]) * meVinK)

                if (words[0:1] == ['pair#']):
                    section = 'pair'
                if (section == 'pair'):
                    if (len(words) == 13):
                        offset = np.array([int(words[6]), int(words[7]), int(words[8])], np.int32)
                        if (self.__valid_offset(offset)):
                            self.bonds.append([int(words[1]), int(words[2]), offset, int(words[12])-1])

                if (len(words) == 0):
                    if (section != 'J_ij'):
                        section = ''
                if (words[0:5] == ['Tc', '(in', 'mean', 'field', 'approximation)']):
                    section = ''

    def __valid_offset(self, offset):
        """
        For removing duplicated bond information (internal use)
        """
        if (offset[0] + offset[1] + offset[2] < 0):
            return False
        elif (offset[0] + offset[1] + offset[2] == 0):
            if (offset[0] + offset[1] < 0):
                return False
            elif (offset[0] + offset[1] == 0):
                if (offset[0] < 0):
                    return False
        return True

def meanfield(unitcell, cutoff):
    """
    Calculate critical temperature by mean-field approximation
    for Ising and Heisenberg models
    """
    # create zJ matrix
    dim = len(unitcell.sites)
    zJ = np.zeros([dim,dim], order='F')
    for source, target, offset, bond_type in unitcell.bonds:
        if abs(unitcell.coupling_constants[bond_type]) > cutoff:
            zJ[source-1][target-1] += unitcell.coupling_constants[bond_type]
            zJ[target-1][source-1] += unitcell.coupling_constants[bond_type]
    # diagonalize zJ matrix
    w, v = np.linalg.eigh(zJ)
    tc = max(w)
    return [tc, tc/3]

if __name__ == '__main__':
    import sys
    if (len(sys.argv) <= 2):
        print('Error: {} path name [cutoff]'.format(sys.argv[0]))
        sys.exit(127)
    path = sys.argv[1]
    name = sys.argv[2]
    cutoff = 0.0
    if (len(sys.argv) > 3):
        cutoff = float(sys.argv[3])
    output_xml = "lattice_" + name + ".xml"
    output_dat = "coupling_" + name + ".dat"
    p = unitcell(path, cutoff)
    print('[[input]]')
    print('path: {}'.format(path))
    print('name: {} {}'.format(output_xml, output_dat))
    print('cutoff: {}'.format(cutoff))
    print('[[unitcell]]')
    print('dimension: {}'.format(p.dimension))
    print('lattice constants: {}'.format(p.lattice_constant))
    print('unit cell basis vectors:')
    for i in  range(p.dimension):
        for j in  range(p.dimension):
            sys.stdout.write('{} '.format(p.cell_bases[i][j]))
        sys.stdout.write('\n')
    print('unit cell volume: {}'.format(p.volume))

    lattices = ElementTree.Element('LATTICES')
    lattice = ElementTree.SubElement(lattices, 'LATTICE')
    lattice.set('name', name)
    lattice.set('dimension', str(p.dimension))
    basis = ElementTree.SubElement(lattice, 'BASIS')
    for j in range(p.dimension):
        vector = ElementTree.SubElement(basis, 'VECTOR')
        el = []
        for i in range(p.dimension):
            el.append(str(p.cell_bases[i][j]))
        vector.text = ' '.join(el)

    unitcell = ElementTree.SubElement(lattices, 'UNITCELL')
    unitcell.set('name', name)
    unitcell.set('dimension', str(p.dimension))
    for s in p.sites:
        vertex = ElementTree.SubElement(unitcell, 'VERTEX')
        vertex.set('type', str(s[1]))
        coordinate = ElementTree.SubElement(vertex, 'COORDINATE')
        el = []
        for i in range(p.dimension):
            el.append(str(s[0][i]))
        coordinate.text = ' '.join(el)
    num_bonds = 0
    for b in p.bonds:
        if (abs(p.coupling_constants[b[3]]) > cutoff):
            edge = ElementTree.SubElement(unitcell, 'EDGE')
            edge.set('type', str(b[3]))
            source = ElementTree.SubElement(edge, 'SOURCE')
            source.set('vertex', str(b[0]))
            target = ElementTree.SubElement(edge, 'TARGET')
            target.set('vertex', str(b[1]))
            el = []
            for i in range(p.dimension):
                el.append(str(b[2][i]))
            target.set('offset', ' '.join(el))
            num_bonds += 1

    print('number of sites: {}'.format(len(p.sites)))
    print('number of bonds: {}'.format(num_bonds))

    root = ElementTree.ElementTree(lattices)
    with open(output_xml, mode='w') as f:
        f.write(minidom.parseString(ElementTree.tostring(lattices, 'utf-8')).toprettyxml(indent="  "))

    with open("coupling_" + name + ".dat", mode='w') as f:
        for i in range(len(p.coupling_constants)):
            if (abs(p.coupling_constants[i]) > cutoff):
                print('  type {}: J = {}'.format(i, p.coupling_constants[i]))
                f.write('{} {}\n'.format(i, p.coupling_constants[i]))

    tc0, tc1 = meanfield(p, cutoff)
    print('Tc (Ising) =', tc0)
    print('Tc (Heisenberg) =', tc1)
