# -*- coding: utf-8 -*-
"""
Ejemplo 16.2 del texto:
Análisis Estructural
R.C. Hibbeler.
Octava Edición.
Pearson.

Este ejemplo utiliza vigas con deformación axial.
"""

import numpy as np
import estructuras as est


def readin():
    nodes    = np.loadtxt('files/H_nodes.txt', ndmin=2)
    mats     = np.loadtxt('files/H_mater.txt', ndmin=2)
    elements = np.loadtxt('files/H_eles.txt' , ndmin=2)
    loads    = np.loadtxt('files/H_loads.txt', ndmin=2)
    return nodes, mats, elements, loads


nodes, mats, elements, loads = readin()
DME_mat, IBC, neq = est.DME(nodes, elements)
KG = est.assembly(elements, mats, nodes, neq, DME_mat)
RHSG = est.loadasem(loads, IBC, neq, 1)
UG = np.linalg.solve(KG, RHSG)
print(UG)