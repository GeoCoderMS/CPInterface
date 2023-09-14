#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import math
import sys
from glob import glob

print(sys.argv)
len_fault_x=float(sys.argv[1])
len_fault_y=float(sys.argv[2])
fault_x = int(sys.argv[3])
fault_y = int(sys.argv[4])
fault_num = fault_x *fault_y
fault_mid_x=float(sys.argv[5])
fault_mid_y=float(sys.argv[6])
fault_mid_z=float(sys.argv[7])

domain_x = float(sys.argv[8])
domain_y = float(sys.argv[9])
domain_z = float(sys.argv[10])

dip=int(sys.argv[11])
strike=float(sys.argv[12])
anomaly = str(sys.argv[13])

home=os.getcwd()


def rotate_point(x_ori,y_ori,cx,cy,angle):
    """Rotate a point around another point by a given angle."""
    angle_rad = math.radians(-angle)
    # Translate point to origin
    x_ori -= cx
    y_ori -= cy
    # Rotate point
    x_new = x_ori * math.cos(angle_rad) - y_ori * math.sin(angle_rad)
    y_new = x_ori * math.sin(angle_rad) + y_ori * math.cos(angle_rad)
    # Translate point back
    x_new += cx
    y_new += cy
    return x_new, y_new
    
    
f=open('mesh.mphtxt','r')
output_file="pylith.mesh"

domain_x = domain_x/2
domain_y = domain_y/2
lines=f.readlines()
for i in range(len(lines)):
    if 'coordinates' in lines[i]:
        #print(i)
        coor=i+1
        #  +1
    elif '# number of element types' in lines[i]:
        #print(i)
        coor_end=i-1
        #print(coor_end)
        
    elif '4 # number of nodes per element' in lines[i]:
        #print(i)
        mesh=i+3
        #  +3
    elif '# Geometric entity indices' in lines[i]:
        #  +3
        meshend=i+1
        #print(meshend)
f.close()

mkmk = glob('mesh_temp')
if mkmk == []:
    os.mkdir('mesh_temp')
os.chdir('mesh_temp')

f=open('coor.txt','w')
a = 0
for i in range(coor,coor_end):
    f.write('%i %s' % (a, lines[i]))
    a += 1
coor_len = a
print(coor_len)
f. close()

f=open('cell.txt','w')
a = 0 
for i in range(mesh,meshend-3):
    f.write('%i %s' % (a, lines[i]))
    a += 1
cell_len = a
print(cell_len)
f.close()

f=open('id.txt','w')
a = 0
for i in range(meshend,meshend+cell_len):
    f.write('%i %s' % (a, lines[i]))
    a += 1 
id_len = a 
print(id_len)
f.close()



#####################################parameter
from math import sin,tan,pi,cos

sol=1e-5




######################################################################################################################################

theta=dip*pi/180

sol=1e-3
fault_number = 1
A=np.loadtxt('coor.txt')
fault_len = []
fault_edge_len = []


cenx=fault_mid_x
ceny=fault_mid_y

for yy in range(0,fault_y,1):
    for xx in range(0,fault_x,1):
        
        f=open('fault%i.txt'%fault_number,'w')
        a=0
        for i in range(len(A)):
            x=A[i,1];
            y=A[i,2];
            z=A[i,3];
            x,y=rotate_point(x,y,cenx,ceny,strike)
            if -len_fault_x/2+((len_fault_x/fault_x)*xx)-sol+fault_mid_x < x <-len_fault_x/2+(len_fault_x/fault_x)*(xx+1)+sol+fault_mid_x and -cos(theta)*len_fault_y/2+cos(theta)*(len_fault_y/fault_y)*(yy)-sol+fault_mid_y<y<-cos(theta)*len_fault_y/2+cos(theta)*(len_fault_y/fault_y)*(yy+1)+sol+fault_mid_y and abs(z+tan(theta)*(y-fault_mid_y)+fault_mid_z)<sol:
                f.write('%i\n'% i)
                a+=1
        fault_len.append(a)
        f.close()

        f=open('fault_edge%i.txt'%fault_number,'w')
        a=0
        for i in range(len(A)):
            x=A[i,1];
            y=A[i,2];
            z=A[i,3];
            x,y=rotate_point(x,y,cenx,ceny,strike)
            if -len_fault_x/2+((len_fault_x/fault_x)*xx)-sol+fault_mid_x < x <-len_fault_x/2+(len_fault_x/fault_x)*(xx+1)+sol+fault_mid_x and -cos(theta)*len_fault_y/2+cos(theta)*(len_fault_y/fault_y)*(yy)-sol+fault_mid_y<y<-cos(theta)*len_fault_y/2+cos(theta)*(len_fault_y/fault_y)*(yy+1)+sol+fault_mid_y and abs(z+tan(theta)*(y-fault_mid_y)+fault_mid_z)<sol:
                
                if  -cos(theta)*len_fault_y/2+cos(theta)*(len_fault_y/fault_y)*(yy)+sol+fault_mid_y>y or y>-cos(theta)*len_fault_y/2+cos(theta)*(len_fault_y/fault_y)*(yy+1)-sol+fault_mid_y or x<(-len_fault_x/2+((len_fault_x/fault_x)*xx))+sol+fault_mid_x  or x>(-len_fault_x/2+((len_fault_x/fault_x)*(xx+1)))-sol+fault_mid_x :
        
                    f.write('%i\n'% i)
                    a+=1
        fault_edge_len.append(a)
        f.close()
        fault_number +=1
#####################################################################################################################################
#####################################################################################################################################



###########################----left----#######################################
sol=1e-5
f=open('left.txt','w')
A=np.loadtxt('coor.txt')
a=0
for i in range(len(A)):
    x=A[i,1];
    y=A[i,2];
    z=A[i,3];
    if x<-domain_x+sol and z>-domain_z+sol or x>domain_x-sol and z>-domain_z+sol :
        f.write('%i \n'% i)
        a+=1
left_len=a
f.close()

###########################----right----#######################################
sol=1e-6
f=open('right.txt','w')
A=np.loadtxt('coor.txt')
a=0
for i in range(len(A)):
    x=A[i,1];
    y=A[i,2];
    z=A[i,3];
    if y>domain_y-sol and z>-domain_z+sol and -domain_x+sol<x<domain_x-sol or y<-domain_y+sol and z>-domain_z+sol and -domain_x+sol<x<domain_x-sol :
        f.write('%i \n'% i)
        a+=1
right_len=a
f.close()



###########################----bottom----#######################################
sol=1e-6
f=open('bottom.txt','w')
A=np.loadtxt('coor.txt')
a=0
for i in range(len(A)):
    x=A[i,1];
    y=A[i,2];
    z=A[i,3];
    if  z<-domain_z+sol  :
        f.write('%i \n'% i)
        a+=1
bottom_len=a
f.close()


###############################################################################################################
################################----test.mesh----##############################################################
###############################################################################################################


os.chdir(home)

f=open('pylith.mesh','w')

f.write('mesh = { \n   dimension = 3\n   use-index-zero = true \n   vertices = {\n    dimension = 3\n    count = %i\n    coordinates = {\n'% coor_len)

C=open(home + '/mesh_temp/coor.txt','r')
lines=C.readlines()
for i in range(len(lines)):
    f.write(lines[i])
C.close()

f.write('    }\n  }\n  cells = {\n    count = %i\n    num-corners = 4\n    simplices = {\n'%cell_len)

M=open(home + '/mesh_temp/cell.txt','r')
lines=M.readlines()
for i in range(len(lines)):
    f.write(lines[i])
M.close()

f.write('    }\n    material-ids = {\n')


I=open(home + '/mesh_temp/id.txt','r')
lines=I.readlines()
for i in range(len(lines)):
    f.write(lines[i])
I.close()





###########################################################fault######################################################################
for ia in range(1,fault_number,1):
    

    f.write('    }\n  }\n  group = {\n    name = fault%i\n    type = vertices\n    count = %i\n    indices = {\n'%(ia,fault_len[ia-1]))
    
    F=open(home + '/mesh_temp/fault%i.txt'%ia,'r')
    lines=F.readlines()
    for i in range(len(lines)):
        f.write(lines[i])
    F.close()

    f.write('    }\n  }\n  group = {\n    name = fault_edge%i\n    type = vertices\n    count = %i\n    indices = {\n'%(ia,fault_edge_len[ia-1]))

    ###########################----fault_edge----#######################################
    FE=open(home + '/mesh_temp/fault_edge%i.txt'%ia,'r')
    lines=FE.readlines()
    for i in range(len(lines)):
        f.write(lines[i])
    FE.close()
    


###########################################################fault######################################################################





f.write('    }\n  }\n  group = {\n    name = bottom\n    type = vertices\n    count = %i\n    indices = {\n'%bottom_len)

Fyz=open(home + '/mesh_temp/bottom.txt','r')
lines=Fyz.readlines()
for i in range(len(lines)):
    f.write(lines[i])
Fyz.close()

f.write('    }\n  }\n  group = {\n    name = left\n    type = vertices\n    count = %i\n    indices = {\n'%left_len)

N=open(home + '/mesh_temp/left.txt','r')
lines=N.readlines()
for i in range(len(lines)):
    f.write(lines[i])
N.close()


f.write('    }\n  }\n  group = {\n    name = right\n    type = vertices\n    count = %i\n    indices = {\n'%right_len)

N=open(home + '/mesh_temp/right.txt','r')
lines=N.readlines()
for i in range(len(lines)):
    f.write(lines[i])
N.close()

f.write('    }\n  }\n}\n')

f.close()





#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
for i in range(1,fault_num+1):
    CFG = open("fault_%i.cfg"%i,'w')
    CFG.write('''[pylithapp]
[pylithapp.journal.info]
quadrature2d = 1
faultcohesivekin = 1


[pylithapp.timedependent]
formulation = pylith.problems.Implicit
bc = [bottom,left,right]
interfaces = [fault]


[pylithapp.timedependent.bc.bottom]

bc_dof = [0,1,2]

label = bottom

db_initial.label = Dirichlet BC


[pylithapp.timedependent.bc.left]
bc_dof = [0]

label = left

db_initial.label = Dirichlet BC


[pylithapp.timedependent.bc.right]
bc_dof = [1]

label = right

db_initial.label = Dirichlet BC


###############################fault#########################################
[pylithapp.timedependent.interfaces]
[pylithapp.timedependent.interfaces.fault]
id=101
label = fault%i
edge = fault_edge%i
quadrature.cell.dimension = 2
[pylithapp.timedependent.interfaces.fault.eq_srcs.rupture.slip_function]
slip.label = Final slip
slip.iohandler.filename = dislocation_slip.spatialdb
slip_time.label = Slip time
slip_time.iohandler.filename = dislocation_sliptime.spatialdb


# ----------------------------------------------------------------------
# output
# ----------------------------------------------------------------------
# Give basename for VTK output of solution over domain.
[pylithapp.problem.formulation.output.output.writer]
filename = output/dislocation%s.vtk
time_constant = 2*s

#[pylithapp.problem.formulation.output.output]
#vertex_data_fields = [displacement,velocity]


[pylithapp.timedependent.materials.material.output.writer]
filename = output_fault/domain.vtk
time_constant = 2*s'''%(i,i,i))

    if anomaly=='on':
        CFG.write('''
[pylithapp.timedependent.materials.anomaly.output.writer]
filename = output_fault/anomaly.vtk
time_constant = 2*s
''')

    CFG.write('''
[pylithapp.timedependent.interfaces.fault.output.writer]
filename = output_fault/fault%s.vtk
time_constant = 2*s



'''%(i))
    CFG.close()
    print(i)
#########################################################################################
#########################################################################################

F=open("pylithapp.cfg",'w')
F.write('''
[pylithapp]

# ----------------------------------------------------------------------
# journal
# ----------------------------------------------------------------------
[pylithapp.journal.info]
timedependent = 1
explicit = 1
implicit = 1
petsc = 1
solverlinear = 1
meshioascii = 1
homogeneous = 1
implicitelasticity = 1
quadrature2d = 1
fiatsimplex = 1

# ----------------------------------------------------------------------
# mesh_generator
# ----------------------------------------------------------------------
[pylithapp.mesh_generator]
debug = 1

[pylithapp.mesh_generator.reader]
filename = pylith.mesh
coordsys.space_dim = 3

# ----------------------------------------------------------------------
# problem
# ----------------------------------------------------------------------
[pylithapp.problem.formulation]
[pylithapp.timedependent]
dimension = 3
normalizer.length_scale = 1.0*m

[pylithapp.timedependent.formulation.time_step]
total_time = 2*s
dt = 2*s

# ----------------------------------------------------------------------
# materials
# ----------------------------------------------------------------------
[pylithapp.timedependent]
materials = [material''')

if anomaly=='on':
    F.write(',anomaly')
    
F.write(''']

[pylithapp.timedependent.materials]
material = pylith.materials.ElasticIsotropic3D
''')

if anomaly=='on':
    F.write("anomaly = pylith.materials.ElasticIsotropic3D\n")

F.write('''
[pylithapp.timedependent.materials.material]
label = domain

id = 1

db_properties.label = Elastic properties
db_properties.iohandler.filename = matprops.spatialdb

quadrature.cell.dimension = 3
''')

if anomaly =='on':
    F.write('''
[pylithapp.timedependent.materials.anomaly]
label = anomaly

id = 2

db_properties.label = Elastic properties
db_properties.iohandler.filename = anomaly.spatialdb

quadrature.cell.dimension = 3
''')

F.write('''
# ----------------------------------------------------------------------
# PETSc
# ----------------------------------------------------------------------
[pylithapp.problem.formulation]
split_fields = True
use_custom_constraint_pc = True
matrix_type = aij


[pylithapp.petsc]
#pc_type = jacobi
fs_pc_type = fieldsplit
fs_pc_use_amat = true
fs_pc_fieldsplit_type = multiplicative
fs_fieldsplit_displacement_pc_type = ml
fs_fieldsplit_lagrange_multiplier_pc_type = jacobi
fs_fieldsplit_displacement_ksp_type = preonly
fs_fieldsplit_lagrange_multiplier_ksp_type = preonly
''')

F.close()

#########################################################################################
#########################################################################################
#########################################################################################

    
    
    
    

