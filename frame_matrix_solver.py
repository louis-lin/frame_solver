# -*- coding: utf-8 -*-
"""
Created on Thu Dec 31 15:14:27 2020

@author: Louis Lin
"""

import numpy as np

# %%Collections

# ELement collections
class element_collection:
    def __init__(self):
        pass    
    def add_element(self, name,section,node_1,node_2,rez_1 = None ,rez_2=None):
        new_element = element(name,section,node_1,node_2,rez_1,rez_2)
        setattr(self, new_element.name, new_element)

# Section Collection
class section_collection:
    def __init__(self):
        pass
    def add_section(self, name, E, I, A = 1):
        new_section = section(name, E, I , A)
        setattr(self, new_section.name, new_section)
        
# Node Collection
class nodal_collection:
    def __init__(self):
        pass   
    def add_node(self, name, x, y):
        new_node = node(self, name, x, y)
        setattr(self, new_node.name, new_node)       
    def set_property(self, name,**kwargs):
        node = getattr(self, name)
        [setattr(node, str(attr), value) for attr, value in kwargs.items()]  
        # NEED TO ADD FUNCTIONALITY WHERE LOAD IS NOT REPLACE WHEN NEW ONES ARE ADDED
        
# %% Objects
# Define section properties used in the system
class section():
    def __init__(self, name, E, I, A):
        self.name, self.E, self.I, self.A= (name, E, I, A)   
        
# Define nodes in the system
class node:
    def __init__(self,system, name,x,y):
        self.name, self.x, self.y = name, x, y
        self.id = [(len(vars(system)))*[3,3,3][i] + [0,1,2][i] for i in range(3)]
        self.fx , self.fy, self.mz = 0, 0, 0
        self.ux , self.uy, self.uz = 1, 1, 1    
    def get_nodal_forces(self):
        return [self.fx,self.fy,self.mz]       
    def get_free_dof(self):
        return [self.id[i] for i in range(3) if [self.ux, self.uy, self.uz][i]]    
        
# Defining elements in the system
class element:
    def __init__(self,name,section,node_1,node_2, rez_1 ,rez_2 ):
        # Create element with the properties passed in
        self.name = name
        self.node_i , self.node_j = node_1, node_2
        self.E, self.I, self.A = section.E ,section.I, section.A
        
        # Define some properties
        self.u_rigid = self.u_global = self.u_local = np.zeros((6,1))
        self.f_rigid = self.f_global = self.f_local = np.zeros((6,1))
        self.ff_rigid = self.ff_global = self.ff_local = np.zeros((6,1))
        
        # Doesn't always have to pass in a rigid end zone node
        if rez_1 == None and rez_2 == None:
            self.rez_i, self.rez_j = node_1, node_2
        else:
            self.rez_i, self.rez_j = rez_1, rez_2
            
        self.id = self.rez_i.id + self.rez_j.id
        dx = self.node_j.x - self.node_i.x
        dy = self.node_j.y - self.node_i.y
        self.L = np.sqrt(dx**2 + dy**2)
        
        # Define transformation matrices
        self.rotation_transformation(dx,dy,self.L)
        self.rigid_end_transformation()

        self.k_local = self.local_stiffness(self.A,self.E,self.I,self.L)
        self.k_global = np.dot(np.dot(np.transpose(self.g_rot),self.k_local),self.g_rot)
        self.k_rigid = np.dot(np.dot(np.transpose(self.g_rez),self.k_global),self.g_rez)    
    
    # Allows for the definition of the end forces
    def set_fixed_end_forces(self, fx1=0 ,fy1=0, mz1=0, fx2= 0, fy2=0, mz2=0):
        self.ff_local = np.array([[fx1, fy1, mz1, fx2, fy2, mz2]]).transpose()

    def add_fixed_end_forces(self, fx1=0 ,fy1=0, mz1=0, fx2= 0, fy2=0, mz2=0):
        self.ff_local = self.ff_local + np.array([[fx1, fy1, mz1, fx2, fy2, mz2]]).transpose() 
        
    # Defines rotation transformation matrix   
    def rotation_transformation(self,delta_x,delta_y,length):
        cos_phi = np.divide(delta_x ,length, out = np.zeros(1), where= delta_x!=0) # Resolves division by 0, using numpy
        sin_phi = np.divide(delta_y ,length, out = np.zeros(1), where= delta_y!=0)
        ROT = np.array(
                      [[cos_phi[0], sin_phi[0], 0],
                      [-sin_phi[0], cos_phi[0], 0],
                      [0,           0,          1]])
        self.g_rot = np.block([[ROT, np.zeros(ROT.shape)],
                               [np.zeros(ROT.shape),ROT]]) 
        
    # Defines rigid end zone transformation matrix    
    def rigid_end_transformation(self):
        dx_1 = self.node_i.x - self.rez_i.x
        dy_1 = self.node_i.y - self.rez_i.y
        dx_2 = self.node_j.x - self.rez_j.x
        dy_2 = self.node_j.y - self.rez_j.y
        REZ1 = np.array(
                      [[1,  0, -dy_1],
                      [0, 1, dx_1],
                      [0, 0, 1]])   
        REZ2 = np.array(
                      [[1,  0, -dy_2],
                      [0, 1, dx_2],
                      [0, 0, 1]])
        self.g_rez = np.block([[REZ1, np.zeros(REZ1.shape)],
                               [np.zeros(REZ2.shape),REZ2]])
    
    def local_stiffness(self,A,E,I,L):        
        return np.array([[E*A/L, 0,           0,          -E*A/L, 0,           0 ],
                            [0,     12*E*I/L**3, 6*E*I/L**2, 0 ,     -12*E*I/L**3, 6*E*I/L**2],
                            [0,     6*E*I/L**2,  4*E*I/L,    0,      -6*E*I/L**2,  2*E*I/L],
                            [-E*A/L,0,           0,          E*A/L,  0,            0 ],
                            [0,     -12*E*I/L**3,-6*E*I/L**2,0 ,     12*E*I/L**3, -6*E*I/L**2],
                            [0,     6*E*I/L**2,  2*E*I/L,    0,      -6*E*I/L**2,  4*E*I/L]
                            ])
         
def matprint(mat, fmt="g"):
    col_maxes = [max([len(("{:"+fmt+"}").format(x)) for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="  ")
        print("")
# %% Solver
class solve:
    def __init__(self, element_collection, node_collection): 
        self.elements = element_collection
        self.nodes = node_collection
        
        self.no_dof = self.count_dof()
        self.free_dof = self.get_free_dof_id()
        self.support_dof = np.setdiff1d(np.arange(self.no_dof),self.free_dof)
        
        self.u_struct = np.zeros((self.no_dof,1))        
        self.k_struct = self.get_k_struct()
        self.f_struct = self.get_f_struct()       
        
        self.f_free = self.get_partition(self.f_struct, self.free_dof)
        self.k_free_free = self.get_partition(self.k_struct, self.free_dof,self.free_dof)
        self.u_free = np.matmul(np.linalg.inv(self.k_free_free),self.f_free)
        
        self.k_support_free = self.get_partition(self.k_struct, self.support_dof,self.free_dof)
        self.f_support = np.dot(self.k_support_free, self.u_free)
        self.set_from_partition(self.u_struct, self.u_free, self.free_dof)
                
        self.set_element_displacement()
        self.calc_element_forces()
        
        
    def count_dof(self):
        no_dof = max(max([getattr(self.nodes,n).id for n in self.nodes.__dict__.keys()]))+1
        return no_dof
    
    def get_free_dof_id(self):
        free_dofs = []
        for nod in self.nodes.__dict__.keys():
            dof = getattr(self.nodes,nod).get_free_dof()
            free_dofs.append(dof) 
        free_dofs = [dof for sets in free_dofs for dof in sets]
        return free_dofs
    
    def get_partition(self, matrix, rows, cols= 0):
        if cols == 0:
          return [matrix[id_row] for id_row in rows]
        else: 
            partition = np.zeros((len(rows), len(cols)))
            for row, id_row in enumerate(rows):
                for col, id_col in enumerate(cols):
                    partition[row, col] = matrix[id_row, id_col]
            return partition
    
    def set_from_partition(self, matrix, partition, rows, cols = 0):
        if cols == 0:
            for row, id_row in enumerate(rows):
                matrix[id_row] = partition[row]
        else:
            for row, id_row in enumerate(rows):
                for col, id_col in enumerate(cols):
                    matrix[id_row, id_col] += partition[row, col]
    
    def get_k_struct(self):
        k_struct = np.zeros((self.no_dof ,self.no_dof))
        for beam in self.elements.__dict__.keys():
            k = getattr(self.elements,beam).k_rigid
            ids = getattr(self.elements,beam).id
            self.set_from_partition(k_struct,k, rows=ids, cols=ids)
        return k_struct
    
    def get_f_struct(self):
        structural_force = np.zeros((self.no_dof,1))
        for node_name in self.nodes.__dict__.keys():        
            n = getattr(self.nodes,node_name)
            f = n.get_nodal_forces()
            i = 0
            for idx in n.id:
                structural_force[idx] += f[i]
                i += 1
        for beam_name in self.elements.__dict__.keys():
            element = getattr(self.elements, beam_name)
            element.ff_global = np.dot(element.g_rot.transpose(),element.ff_local)
            element.ff_rigid = np.dot(element.g_rez.transpose(),element.ff_global)
            i = 0
            for idx in element.id:
                structural_force[idx] += - element.ff_rigid[i] 
                i += 1
        return structural_force
        
    def set_element_displacement(self):
        for element_name in self.elements.__dict__.keys():
            element = getattr(self.elements, element_name)
            ids = element.id
            for i in range(len(ids)):
                element.u_rigid[i] = self.u_struct[ids[i]]
            element.u_global = np.dot(element.g_rez,element.u_rigid)
            element.u_local = np.dot(element.g_rot,element.u_global)
        
    def calc_element_forces(self):
        for element_name in self.elements.__dict__.keys():
            element = getattr(self.elements, element_name)            
            for system in ['_rigid', '_global', '_local']:
                k = getattr(element, "k"+ system)
                u = getattr(element, "u"+ system)
                ff = getattr(element, "ff" + system)
                setattr(element, "f"+system, np.dot(k,u) + ff)
                