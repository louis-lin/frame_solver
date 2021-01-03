# -*- coding: utf-8 -*-
"""
Created on Fri Jan  1 01:04:28 2021

@author: Louis Lin
"""

import frame_matrix_solver as solver
from math import sqrt

node = solver.nodal_collection()
element = solver.element_collection()
section = solver.section_collection()

E = 2*10**11
section.add_section("ab",E, 2*10**-4, 6*10**-3)   
section.add_section("bc",E, 5*10**-5, 4*10**-3)   

node.add_node("A",0,0)
node.add_node("B",8,0)
node.add_node("C",13,0)


element.add_element("Beam_AB",section.ab,node.A,node.B)
element.add_element("Beam_BC",section.bc,node.B,node.C)

P = 5/sqrt(2) * 10**3

node.set_property("C",fx = P, fy = -P)

node.set_property("A",ux = 0, uy = 0, uz = 0)
node.set_property("B", uy = 0)
# node.set_property("C",ux = 0, uy = 0, uz = 0)

solution = solver.solve(element,node)