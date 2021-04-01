"""
Modulo de funciones para realizar analisis de estructuras
aporticadas por el metodo de rigidez.
---------------------------------------------------------
Para uso en la materia
IC0285 Modelacion Computacional
del pregrado de Ingenieria Civil de la
UNIVERSIDAD EAFIT.
---------------------------------------------------------
Se permiten elementos tipo viga con y sin deformacion
axial.


@author: Juan David Gomez
"""

import numpy as np


def readin():
    """Lee los archivos del modelo y los carga en arreglos

    Retorna
    -------
    nodes:    nd array
      Arreglo con la informacion de los nodos del modelo (nnodes , 6)
    elements: nd array
      Arreglo con la informacion de los elementos del modleo (numel , 5)
    mats: nd array
      Arreglo con los perfiles de material.
    loads: nd array
      Arreglo con la informacion de las cargas nodales del modelo (nloads , 3)
    
    """
    nodes = np.loadtxt('Vnodes.txt', ndmin=2)
    mats = np.loadtxt('Vmater.txt', ndmin=2)
    elements = np.loadtxt('Veles.txt', ndmin=2)
    loads = np.loadtxt('Vloads.txt', ndmin=2)
    return nodes, mats, elements, loads

def DME(nodes, elements):
    """Calcula la matriz ensambladora de ecuaciones DME()
    Parametros
    ----------
    nodes:    nd array
      Arreglo con la informacion de los nodos del modelo (nnodes , 6)
    elements: nd array
      Arreglo con la informacion de los elementos del modleo (numel , 5)
    Retorna
    -------
    DME_mat : nd array
     Matriz indicadora de ecuaciones (numel , 6)
    IBC     : nd array
     Areglo de condiciones de ecuaciones nodales
    neq     : integer
     Numero total de ecuaciones del modelo.
    
    """

    nels = elements.shape[0]
    IELCON = np.zeros([nels, 2], dtype=np.integer)
    DME_mat = np.zeros([nels, 6], dtype=np.integer)
    neq, IBC = eqcounter(nodes)
    nnodes = 2
    for i in range(nels):
        for j in range(nnodes):
            IELCON[i, j] = elements[i, j+3]
            kk = IELCON[i, j]
            for l in range(3):
                DME_mat[i, 3*j+l] = IBC[kk, l]
    return DME_mat, IBC, neq


def eqcounter(nodes):
    """Identifica ecuaciones activas y las cuenta
    
    Parametros
    ----------
    nodes:    nd array
      Arreglo con la informacion de los nodos del modelo (nnodes , 6)
  
    Retorna
    -------
    neq     : integer
     Numero total de ecuaciones del modelo.
    IBC     : nd array
     Areglo de condiciones de ecuaciones nodales
    
    """
    nnodes = nodes.shape[0]
    IBC = np.zeros([nnodes, 3], dtype=np.integer)
    for i in range(nnodes):
        for k in range(3):
            IBC[i, k] = int(nodes[i, k+3])
    neq = 0
    for i in range(nnodes):
        for j in range(3):
            if IBC[i, j] == 0:
                IBC[i, j] = neq
                neq = neq + 1
    return neq, IBC

def assembly(elements, mats, nodes, neq, DME_mat):
    
    """Calcula la matriz de rigidez global de la estructura KG(neq , neq) para
    ensamblajes de elementos tipo viga con (elemento tipo 0) y sin deformacion
    axial (elemento tipo 1).
    
    Parametros
    ----------
    elements: nd array
      Arreglo con la informacion de los elementos del modleo (numel , 5)
    mats: nd array
      Arreglo con los perfiles de material
    nodes:    nd array
      Arreglo con la informacion de los nodos del modelo (nnodes , 6)
    neq     : integer
     Numero total de ecuaciones del modelo.
    DME_mat : nd array
     Matriz indicadora de ecuaciones (numel , 6)

    Retorna
    -------
    KG     : nd array
     Matriz de rigidez global de la estructura
    
    """
    
    IELCON = np.zeros([2], dtype=np.integer)
    KG = np.zeros((neq, neq))
    nels = elements.shape[0]
    nnodes = 2
    ndof = 6
    for el in range(nels):
        elcoor = np.zeros([nnodes, 2])
        im = np.int(elements[el, 2])
        iet = np.int(elements[el, 1])
        par0 = mats[im , 0]   # Iz
        par1 = mats[im , 1]   # Emod
        par2 = mats[im , 2]   # A
        for j in range(nnodes):
            IELCON[j] = elements[el, j+3]
            elcoor[j, 0] = nodes[IELCON[j], 1]
            elcoor[j, 1] = nodes[IELCON[j], 2]
        if iet == 0:
            kloc = uelbeam2D(elcoor, par0, par1 , par2)
        elif iet == 1:
            kloc = uelbeam2DU(elcoor, par0 , par1)
        dme = DME_mat[el, :ndof]
        for row in range(ndof):
            glob_row = dme[row]
            if glob_row != -1:
                for col in range(ndof):
                    glob_col = dme[col]
                    if glob_col != -1:
                        KG[glob_row, glob_col] = KG[glob_row, glob_col] +\
                            kloc[row, col]
    return KG

def local_forces(elements, mats, nodes, neq, DME_mat , UC):
    
    """Calcula las fuerzas en los elementos en coordenadas locales
    
    Parametros
    ----------
    elements: nd array
      Arreglo con la informacion de los elementos del modleo (numel , 5)
    mats: nd array
      Arreglo con los perfiles de material
    nodes:    nd array
      Arreglo con la informacion de los nodos del modelo (nnodes , 6)
    neq     : integer
     Numero total de ecuaciones del modelo.
    DME_mat : nd array
     Matriz indicadora de ecuaciones (numel , 6)
    UC      : nd array
     Vector de desplazmientos nodales incluyendo las restricciones

    Retorna
    -------
    FG     : nd array
     Vector con las fuerzas de cada elemento en coordenadas locales
    
    """
    IELCON = np.zeros([2], dtype=np.integer)
    nels = elements.shape[0]
    nnodes = 2
#
    for el in range(nels):
        iet = np.int(elements[el , 1])
        if iet == 0:
            ndof = 6
            FG = np.zeros((nels, 6))
            ul = np.zeros(6)
            fl = np.zeros(6)
        elif iet == 1:
            ndof = 4
            FG = np.zeros((nels, 4))
            ul = np.zeros(4)
            fl = np.zeros(4)   
#
    for el in range(nels):
#
        iet = np.int(elements[el , 1])         
#
        elcoor = np.zeros([nnodes, 2])
        im = np.int(elements[el , 2])
        par0 = mats[im , 0]   # Iz
        par1 = mats[im , 1]   # Emod
        par2 = mats[im , 2]   # A
        for j in range(nnodes):
            IELCON[j] = elements[el , j+3]
            elcoor[j, 0] = nodes[IELCON[j] , 1]
            elcoor[j, 1] = nodes[IELCON[j] , 2]        
        for j in range(ndof):
            ig = DME_mat[el, j]
            ul[j] = UC[ig]        
        if iet == 0:            
            fl = reac_beam2D(elcoor , par0, par1 , par2 , ul)
        elif iet == 1:            
            fl = reac_beam2DU(elcoor , par0, par1 , ul)
        FG[el , :] = fl[:]
                
    return FG
#
def global_forces(elements, mats, nodes, neq, DME_mat , UC):
    """Calcula las fuerzas en los elementos en coordenadas globales
    
    Parametros
    ----------
    elements: nd array
      Arreglo con la informacion de los elementos del modleo (numel , 5)
    mats: nd array
      Arreglo con los perfiles de material
    nodes:    nd array
      Arreglo con la informacion de los nodos del modelo (nnodes , 6)
    neq     : integer
     Numero total de ecuaciones del modelo.
    DME_mat : nd array
     Matriz indicadora de ecuaciones (numel , 6)
    UC      : nd array
     Vector de desplazmientos nodales incluyendo las restricciones

    Retorna
    -------
    FG     : nd array
     Vector con las fuerzas de cada elemento en coordenadas globales
    
    """
    IELCON = np.zeros([2], dtype=np.integer)
    nels = elements.shape[0]
    nnodes = 2
#
    for el in range(nels):
        iet = np.int(elements[el , 1])
        if iet == 0:
            ndof = 6
            FG = np.zeros((nels, 6))
            ul = np.zeros(6)
            fl = np.zeros(6)
        elif iet == 1:
            ndof = 4
            FG = np.zeros((nels, 4))
            ul = np.zeros(4)
            fl = np.zeros(4)   
#
    for el in range(nels):
#
        iet = np.int(elements[el , 1])         
#
        elcoor = np.zeros([nnodes, 2])
        im = np.int(elements[el , 2])
        par0 = mats[im , 0]   # Iz
        par1 = mats[im , 1]   # Emod
        par2 = mats[im , 2]   # A
        for j in range(nnodes):
            IELCON[j] = elements[el , j+3]
            elcoor[j, 0] = nodes[IELCON[j] , 1]
            elcoor[j, 1] = nodes[IELCON[j] , 2]        
        for j in range(ndof):
            ig = DME_mat[el, j]
            ul[j] = UC[ig]        
        if iet == 0:            
            fl = reac_beam2D_global(elcoor , par0, par1 , par2 , ul)
        elif iet == 1:            
            fl = reac_beam2DU_global(elcoor , par0, par1 , ul)
        FG[el , :] = fl[:]
                
    return FG
#
#
def uelbeam2DU(coord, I, Emod):
    """Elemento viga (2D) sin deformacion axial

    Parametros
    ----------
    coord : ndarray
      Cordenadas nodales (2, 2).
    I : float
      Momento de inercia de la seccion transversal.
    Emod : float
      Modulo de elasticidad (>0).

    Retorna
    -------
    kl : ndarray
      Matriz de rigidez local para el elemento (4, 4) rotada al sistema
      global de coordenadas.

    """
    vec = coord[1, :] - coord[0, :]
    nx = vec[0]/np.linalg.norm(vec)
    ny = vec[1]/np.linalg.norm(vec)
    L = np.linalg.norm(vec)
    Q = np.array([
        [-ny, nx, 0,  0,  0, 0],
        [0,  0, 1.0,  0,  0, 0],
        [0,  0, 0, -ny, nx, 0],
        [0,  0, 0,  0,  0, 1.0]])
    kl = (I*Emod/(L*L*L)) * np.array([
        [12.0, 6, -12.0, 6*L],
        [6,  4*L*L, -6*L, 2*L*L],
        [-12.0,  -6*L, 12.0, -6*L],
        [6*L,  2*L*L, -6*L, 4*L*L]])
    kG = np.dot(np.dot(Q.T, kl), Q)
    return kG
#
def loadasem(loads, IBC, neq):
    """Ensambla el vector de cargas globales RHSG(neq)

    Parametros
    ----------
    loads: nd array
      Arreglo con la informacion de las cargas nodales del modleo (nloads , 3)
    IBC     : nd array
     Areglo de condiciones de ecuaciones nodales
    neq     : integer
     Numero total de ecuaciones del modelo.


    Retorna
    -------
    RHSG : ndarray
      vector de cargas globales RHSG(neq).

    """
    RHSG = np.zeros([neq])
    nl = loads.shape[0]
    for i in range(nl):
        il = int(loads[i, 0])
        ilx = IBC[il, 0]
        ily = IBC[il, 1]
        ilT = IBC[il, 2]
        if ilx != -1:
            RHSG[ilx] = loads[i, 1]
        if ily != -1:
            RHSG[ily] = loads[i, 2]
        if ilT != -1:
            RHSG[ilT] = loads[i, 3]
    return RHSG
#
def uelbeam2D(coord, Iz , Emod , A):
    """Elemento viga (2D) con deformacion axial

    Parametros
    ----------
    coord : ndarray
      Cordenadas nodales (2, 2).
    Iz : float
      Momento de inercia de la seccion transversal.
    Emod : float
      Modulo de elasticidad (>0).
    A    : float
      Area de la sección transversal de la viga.

    Retorna
    -------
    kl : ndarray
      Matriz de rigidez local para el elemento (6 , 6) rotada al sistema
      global de coordenadas.

    """
    vec    = coord[1, :] - coord[0, :]
    L      = np.linalg.norm(vec)
    Q      = Q_2Dframe(vec)
    #
    kl = (Iz*Emod/(L*L*L)) * np.array([[A*L*L/Iz   , 0.0  , 0.0   , -A*L*L/Iz , 0.0   , 0.0  ],
                                       [ 0.0       , 12.0 , 6*L   , 0.0       , -12.0 , 6*L  ],
                                       [ 0.0       , 6*L  , 4*L*L , 0.0       , -6*L  , 2*L*L],
                                       [-A*L*L/Iz  , 0.0  , 0.0   , A*L*L/Iz  , 0.0   , 0.0  ],
                                       [0.0        ,-12.0 , -6*L  , 0.0       , 12.0  , -6*L ],
                                       [0.0        , 6*L  , 2*L*L , 0.0       , -6*L  , 4*L*L]])
    #
    kG = np.dot(np.dot(Q.T, kl), Q)
    
    return kG
    #
def reac_beam2D(coord, Iz , Emod , A , U):
    """Calcula las fuerzas internas en coordenadas locales
       para una viga con deformacion axial

    Parametros
    ----------
    coord : ndarray
      Cordenadas nodales (2, 2).
    Iz : float
      Momento de inercia de la seccion transversal.
    Emod : float
      Modulo de elasticidad (>0).
    A    : float
      Area de la seccion transversal de la viga.
    U    : nd array
      Vector con desplazamientos nodales para el elemento
      en coordenadas globales

    Retorna
    -------
    fl : ndarray
      Vector de fuerzas internas en coordenadas locales.

    """
    vec    = coord[1, :] - coord[0, :]
    L      = np.linalg.norm(vec)
    Q      = Q_2Dframe(vec)
    ul = np.dot(Q , U)
    #
    kl = (Iz*Emod/(L*L*L)) * np.array([[A*L*L/Iz   , 0.0  , 0.0   , -A*L*L/Iz , 0.0   , 0.0  ],
                                       [ 0.0       , 12.0 , 6*L   , 0.0       , -12.0 , 6*L  ],
                                       [ 0.0       , 6*L  , 4*L*L , 0.0       , -6*L  , 2*L*L],
                                       [-A*L*L/Iz  , 0.0  , 0.0   , A*L*L/Iz  , 0.0   , 0.0  ],
                                       [0.0        ,-12.0 , -6*L  , 0.0       , 12.0  , -6*L ],
                                       [0.0        , 6*L  , 2*L*L , 0.0       , -6*L  , 4*L*L]])
    #
    fl = np.dot(kl , ul)

    
    return fl
#
def reac_beam2DU(coord, I, Emod , U):
    """Calcula las fuerzas internas en coordenadas locales
       para una viga sin deformación axial

    Parametros
    ----------
    coord : ndarray
      Cordenadas nodales (2, 2).
    I  : float
      Momento de inercia de la seccion transversal.
    Emod : float
      Modulo de elasticidad (>0).
    U    : nd array
      Vector con desplazamientos nodales para el elemento
      en coordenadas globales

    Retorna
    -------
    fl : ndarray
      Vector de fuerzas internas en coordenadas locales.

    """
    vec = coord[1, :] - coord[0, :]
    L = np.linalg.norm(vec)
    nx = vec[0]/L
    ny = vec[1]/L
    
    Q = np.array([
        [-ny, nx, 0,  0,  0, 0],
        [0,  0, 1.0,  0,  0, 0],
        [0,  0, 0, -ny, nx, 0],
        [0,  0, 0,  0,  0, 1.0]])
    ul = np.dot(Q , U)
    kl = (I*Emod/(L*L*L)) * np.array([
        [12.0, 6, -12.0, 6*L],
        [6,  4*L*L, -6*L, 2*L*L],
        [-12.0,  -6*L, 12.0, -6*L],
        [6*L,  2*L*L, -6*L, 4*L*L]])
    fl = np.dot(kl , ul)
    return fl

def reac_beam2D_global(coord, Iz , Emod , A , U):
    """Calcula las fuerzas internas en coordenadas globales
       para una viga con deformacion axial

    Parametros
    ----------
    coord : ndarray
      Cordenadas nodales (2, 2).
    Iz : float
      Momento de inercia de la seccion transversal.
    Emod : float
      Modulo de elasticidad (>0).
    A    : float
      Area de la seccion transversal de la viga.
    U    : nd array
      Vector con desplazamientos nodales para el elemento
      en coordenadas globales

    Retorna
    -------
    fl : ndarray
      Vector de fuerzas internas en coordenadas globales.

    """
    vec    = coord[1, :] - coord[0, :]
    L      = np.linalg.norm(vec)
    Q      = Q_2Dframe(vec)
    #
    kl = (Iz*Emod/(L*L*L)) * np.array([[A*L*L/Iz   , 0.0  , 0.0   , -A*L*L/Iz , 0.0   , 0.0  ],
                                       [ 0.0       , 12.0 , 6*L   , 0.0       , -12.0 , 6*L  ],
                                       [ 0.0       , 6*L  , 4*L*L , 0.0       , -6*L  , 2*L*L],
                                       [-A*L*L/Iz  , 0.0  , 0.0   , A*L*L/Iz  , 0.0   , 0.0  ],
                                       [0.0        ,-12.0 , -6*L  , 0.0       , 12.0  , -6*L ],
                                       [0.0        , 6*L  , 2*L*L , 0.0       , -6*L  , 4*L*L]])
    #
    kG = np.dot(np.dot(Q.T, kl), Q)
    fl = np.dot(kG , U)
    
    return fl

def reac_beam2DU_global(coord, I, Emod , U):
    """Calcula las fuerzas internas en coordenadas globales
       para una viga sin deformacion axial

    Parametros
    ----------
    coord : ndarray
      Cordenadas nodales (2, 2).
    I  : float
      Momento de inercia de la seccion transversal.
    Emod : float
      Modulo de elasticidad (>0).
    U    : nd array
      Vector con desplazamientos nodales para el elemento
      en coordenadas globales

    Retorna
    -------
    fl : ndarray
      Vector de fuerzas internas en coordenadas globales.

    """
    vec = coord[1, :] - coord[0, :]
    L = np.linalg.norm(vec)
    nx = vec[0]/L
    ny = vec[1]/L
    
    Q = np.array([
        [-ny, nx, 0,  0,  0, 0],
        [0,  0, 1.0,  0,  0, 0],
        [0,  0, 0, -ny, nx, 0],
        [0,  0, 0,  0,  0, 1.0]])

    kl = (I*Emod/(L*L*L)) * np.array([
        [12.0, 6, -12.0, 6*L],
        [6,  4*L*L, -6*L, 2*L*L],
        [-12.0,  -6*L, 12.0, -6*L],
        [6*L,  2*L*L, -6*L, 4*L*L]])
    kG = np.dot(np.dot(Q.T, kl), Q)
    fl = np.dot(kG , U)
    return fl




def Q_2Dframe(vec):
    """ 2D - frame transformation matrix.
    
    Parameters
    ---------- 
    vec: nd array
      Vector formed by initial and final node coordenates.
    Retorna
    -------
    Q : ndarray
      Matriz de transformación bajo rotacion.

    """
    # 
    L  = np.linalg.norm(vec) 
    nx = vec[0]/L
    ny = vec[1]/L
    #
    Q  = np.array([[ nx ,  ny ,   0     ,  0   ,  0 , 0 ],
                   [-ny ,  nx ,   0     ,  0   ,  0 , 0 ],
                   [  0 ,  0  ,   1     ,  0   ,  0 , 0 ],
                   [  0 ,  0  ,   0     , nx   , ny , 0 ],
                   [  0 ,  0  ,   0     ,-ny   , nx , 0 ],
                   [  0 ,  0  ,   0     ,  0   , 0  , 1 ]])
    return Q

def complete_disp(IBC, nodes, UG):
    """Fills the displacement vectors with imposed and computed values.
    Parameters
    ---------

    IBC : ndarray (int)
        IBC (Indicator of Boundary Conditions) indicates if the
        nodes has any type of boundary conditions applied to it.
    UG : ndarray (float)
        Array with the computed displacements.
    nodes : ndarray (float)
        Array with number and nodes coordinates

    Returns
    -------
    UC : (nnodes, 2) ndarray (float)
      Array with the displacements.

    """
    nnodes = nodes.shape[0]
    UC = np.zeros([nnodes, 3], dtype=np.float)
    for row in range(nnodes):
        for col in range(3):
            cons = IBC[row, col]
            if cons == -1:
                UC[row, col] = 0.0
            else:
                UC[row, col] = UG[cons]

    return UC

def empotramiento(W , l ):
    
    """Calcula las fuerzas de empotramiento
    para una viga de luz L con carga uniformemente
    distribuida de intensidad W

    """
    F = W*l/2.0
    M = W*(l**2)/12.0

    return F , M

