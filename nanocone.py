#!/usr/bin/python
# -*- coding:utf-8 -*-


import numpy as np

#Cubic BN
#   3.57
#   0.0 0.5 0.5
#   0.5 0.0 0.5
#   0.5 0.5 0.0
#1 1
# Direct
#  0.00 0.00 0.00 
#  0.25 0.25 0.25
#The first line is treated as a comment line (you should write down the name of the system). The second line provides a universal scaling factor (lattice constant), which is used to scale all lattice vectors and all atomic coordinates (of this value is negative it is interpreted as the total volume of the cell). On the following three lines the three lattice vectors defining the unit cell of the system are given (first line corresponding to the first lattice vector, second to the second, and third to the third). The sixth line supplies the number of atoms per atomic species (one number for each atomic species). The ordering must be consistent with the POTCAR and the INCAR file. The seventh line switches to selective dynamics (only the first character is relevant and must be S or s). This mode allows to provide extra flags for each atom signaling whether the respective coordinate(s) of this atom will be allowed to change during the ionic relaxation. This setting is useful if only certain shells around a defect or layers near a surface should relax.
#

def get_intersect(s):
            """ 
            Returns the point of intersection of the lines passing through a2,a1 and b2,b1.
            a1: [x, y] a point on the first line
            a2: [x, y] another point on the first line
            b1: [x, y] a point on the second line
            b2: [x, y] another point on the second line
            """
            #s = np.vstack([a1,a2,b1,b2])        # s for stacked
            print(s)
            #h = np.hstack((s, np.ones((4, 1)))) # h for homogeneous
            h =s 
            print(h)
            l1 = np.cross(h[0], h[1])           # get first line
            l2 = np.cross(h[2], h[3])           # get second line
            print(l1)
            print(l2)
            x, y, z = np.cross(l1, l2)          # point of intersection
            if z == 0:                          # lines are parallel
               return (float('inf'), float('inf'))
            return x,y,z

def convert_string_to_vector(string):
    eles = string.split()
    vector=[]
    if len(eles) != 3:
        print("invalid vector: " + string)
        return []
    for e in eles:
        try:
            vector.append(float(e))
        except:
            print("invalid vector: "+ string )
            return []
    return vector


def load_POSCAR(path):
    f = ''
    scaling_factor = 1.0
    lattice_vectors=[]
    atoms = []
    try:
        f = open(path,'r')
        lines = f.readlines()
        if len(lines) < 8 :
            print("POSCAR:"+path+" format error: too short")
            return []

        try:
            scaling_factor = float(lines[1])
        except:
            print("invalid scaling factor:" + lines[1] )
            return []

        for i in range(2,5):
            vector = convert_string_to_vector(lines[i].replace('\n',''))
            if vector == [] :
                print("invalid lattice_vectors")
                return [] 
            lattice_vectors.append(vector)
        # convert to numopy.array
        lattice_vectors = np.array(lattice_vectors)
        atom_type = lines[5].replace('\n','')
        atom_num = int(lines[6].replace('\n',''))
        coordinate_type = lines[7].replace('\n','').strip()
        if coordinate_type != "Direct" :
            print("unsupport coordinate type: "+ coordinate_type)
            return []

        for i in range(8, atom_num+8):
            vector = convert_string_to_vector(lines[i].replace('\n',''))
            if vector !=[] :
                atoms.append(vector)
        # convert to numpy.array
        atoms = np.array(atoms)
        return scaling_factor,lattice_vectors,atom_type,atom_num,atoms

    finally:
        if f != '' :
            f.close()

class POSCAR:
    def __init__(self,path):
     self.scaling_factor = 1.0
     self.lattice_vectors=[]
     self.atoms =[]
     self.coordinate_type=""
     self.atom_type=""
     self.atom_num =0
     self.atoms =[]
     self.scaling_factor, self.lattice_vectors, self.atom_type, self.atom_num, self.atoms = load_POSCAR(path)
     

    def gen_cone(self, path):
        vectors = []
        try:
            f = open(path,'r')
            lines = f.readlines()
            if len(lines) != 4 :
                print("cone parameter:"+path+" format error: too short")
                return []

            for i in range(0, 4):
                vector = convert_string_to_vector(lines[i].replace('\n',''))
                if vector == [] :
                    print("invalid cone parameter")
                    return [] 
                vectors.append(vector)
            lines = np.array(vectors)
            vectors = np.dot(lines,pos.lattice_vectors)
            print(get_intersect(vectors))


        finally:
            if f != '' :
                f.close()
            


if __name__ == "__main__":
    pos = POSCAR("./POSCAR")
    real_pos=np.dot(pos.atoms, pos.lattice_vectors)
    print(real_pos)
    pos.gen_cone("./points")
