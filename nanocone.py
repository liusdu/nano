#!/usr/bin/python
# -*- coding:utf-8 -*-


import numpy as np
import sys

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
     
    def get_vectors(self, path) :
        vectors = []
        try:
            f = open(path,'r')
            lines = f.readlines()
            for i in range(len(lines)):
                vector = convert_string_to_vector(lines[i].replace('\n',''))
                if vector == [] :
                    print("invalid cone parameter")
                    return [] 
                vectors.append(vector)
            lines = np.array(vectors)
            vectors = np.dot(lines, self.lattice_vectors)
            return vectors


        finally:
            if f != '' :
                f.close()

    # all atoms are on x-y plane
    def extend_supercells(self, x, y):
        if x % 2 != 1  or y % 2 != 1 :
            print("invalid  factor ", x,y)
            exit(1)
        after = np.array([])
        for i in range(len(self.atoms)):
            self.atoms[i][2] = 0
        for i in range(0, x):
            for j in range(0, y):
                c_x = i - x/2
                c_y  = j - y/2
                after = np.append(after, self.atoms + np.array([c_x,c_y,0]))
                after = after.reshape(-1,3)
        self.atoms = after

    def gen_cone(self, origin_path, lines_path):
        origin = self.get_vectors(origin_path)
        lines = self.get_vectors(lines_path)
        origin_p = (origin[0] + origin[1])/2
        origin_p[2] = 0
        for i in range(len(lines)) :
            lines[i] = lines[i] - origin_p
            lines[i][2] = 0; 
        line1 =  lines[1] - lines[0]
        line2 =  lines[3] - lines[2]
        self.extend_supercells(3,3)
        cut_ang = np.arccos(np.dot(line1,line2)/(np.linalg.norm(line1)*np.linalg.norm(line2)))
        real_pos=np.dot(pos.atoms, pos.lattice_vectors)
        
        base_x = np.array([1,0,0])
        for i in range(len(real_pos)):
            theta = np.arccos(np.dot(real_pos[i],base_x)/np.linalg.norm(real_pos[i]))
            beta =  theta * np.pi *2/(np.pi*2 - cut_ang)
            x = np.linalg.norm(real_pos[i]) * (np.pi*2 - cut_ang)/(np.pi *2) * np.cos(beta)
            y =  np.linalg.norm(real_pos[i]) * (np.pi*2 - cut_ang)/(np.pi *2)  *np.sin(beta)
            if  real_pos[i][1] < 0 :
                y = y * -1
            z = np.linalg.norm(real_pos[i]) *np.sqrt( 1- np.square((np.pi*2 - cut_ang)/(np.pi *2))) 
            real_pos[i][0] = x
            real_pos[i][1] = y
            real_pos[i][2] = z
        real = (real_pos/100)
        real1 = []
        for i in range(len(real)):
            equal= 0
            for j in range(len(real1)):
                if int(real[i][0]*10000) == int(real1[j][0]*10000) and int(real[i][1]*10000) == int(real1[j][1]*10000)  and int(real[i][2]*10000) == int(real1[j][2]*10000):
                        equal = 1
            if equal == 0 :
                if real[i][2] > 1:
                    continue
                real1.append(real[i])
                print("\t%0.4f\t%0.4f\t%0.4f"%(real[i][0] + 0.5, real[i][1] + 0.5,real[i][2]))
            



if __name__ == "__main__":
    
    if len(sys.argv) != 2 :
        print("useage: nanocone.py POSCAR")

    poscar_file = sys.argv[1]
    pos = POSCAR(poscar_file)

    pos.gen_cone("./origin", "lines")
