#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import os
from math import pi, sqrt

########################### Model Setup ######################################
domain_size = [1000e3,1000e3,1000e3]      # Size of Domain [X, Y, Z]

fault_coor = [0, 0, 50e3]                 # Fault Center Coordinate [X, Y, Z]

fault_hei=100e3                           # Hegiht of Fault
fault_wid=100e3                           # Width of Fault

x_num=4                                   # Number of sub-fault along Strike
y_num=4                                   # Number of sub-fault along Dip

strike_angle=45                           # Strike Angle
dip_angle=30                              # Dip Angle

GPS=np.loadtxt('GPS_station.txt')         # GPS_station location file

############################ Anomaly #########################################
anomaly = 'off'                           # "on" = anomaly / "off" = no anomaly

anomaly_size = [200e3, 200e3, 40e3]       # Size of Anomaly [X, Y, Z]
anomaly_coor = [0, 0, -130e3]             # Anomaly Center Coordinate [X, Y, Z]

############################ Properties ######################################

domain_density = 2700;  anomaly_density = 3000 ;    # Density 
domain_nu      = 0.25;  anomaly_nu      = 0.25 ;    # Nu      (Possion's ratio)
domain_E       = 67e9;  anomaly_E       = 130e9;    # E       (Young's modulus)

##############################################################################

home=os.getcwd()
F=open('green_mesh.java','w')

run_count=0
run=1
F.write('''/*
 * save.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on May 13 2022, 16:36 by COMSOL 5.4.0.346. */
public class green_mesh {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("\ubc14\ud0d5\ud654\uba74");

    model.label("green_mesh.mph");

    model.param().set("domain_X", "%i");
    model.param().set("domain_Y", "%i");
    model.param().set("domain_Z", "%i");
    model.param().set("fault_x", "%i");
    model.param().set("fault_y", "%i");
    model.param().set("fault_z", "%i");
    model.param().set("fault_hei", "%i");
    model.param().set("fault_wid", "%i");
    model.param().set("angle", "%i*pi/180");
    model.param().set("strike", "-90-%f");
    model.param().set("x_num", "%i");
    model.param().set("y_num", "%i");

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom1", 3);

    model.component("comp1").mesh().create("mesh1");

    model.component("comp1").geom("geom1").create("blk1", "Block");
    model.component("comp1").geom("geom1").feature("blk1").set("pos", new String[]{"0", "0", "-domain_Z/2"});
    model.component("comp1").geom("geom1").feature("blk1").set("base", "center");
    model.component("comp1").geom("geom1").feature("blk1")
         .set("size", new String[]{"domain_X", "domain_Y", "domain_Z"});
    model.component("comp1").geom("geom1").create("wp1", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp1").set("planetype", "coordinates");
    model.component("comp1").geom("geom1").feature("wp1")
         .set("genpoints", new String[][]{{"fault_x", "fault_y", "-fault_z"}, {"fault_x+1", "fault_y", "-fault_z"}, {"fault_x", "fault_y+(1/tan(angle))", "-fault_z-1"}});
    model.component("comp1").geom("geom1").feature("wp1").set("unite", true);
'''%(domain_size[0],domain_size[1],domain_size[2],fault_coor[0],fault_coor[1],fault_coor[2],fault_hei,fault_wid,dip_angle,strike_angle,x_num,y_num))
#######################################################################
a=0
for i in range(y_num):
    for j in range(x_num):
        a+=1
        F.write('''
    model.component("comp1").geom("geom1").feature("wp1").geom().create("r%i", "Rectangle");
    model.component("comp1").geom("geom1").feature("wp1").geom().feature("r%i")
         .set("pos", new String[]{"-fault_wid/2+fault_wid/%s*%s+fault_wid/2/%s", "fault_hei/2-fault_hei/%s*%s-fault_hei/2/%s"});
    model.component("comp1").geom("geom1").feature("wp1").geom().feature("r%i").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp1").geom().feature("r%i")
         .set("size", new String[]{"fault_wid/%s", "fault_hei/%s"});
        '''%(a,a,x_num,j,x_num,y_num,i,y_num,a,a,x_num,y_num))
        run_count+=1
        if run_count==60:
            run_count=0
            run+=1
            F.write('''    return model;
  }

  public static Model run%s(Model model) {
            '''%run)

######################################################################
point_num=1

for i in range(len(GPS)):    
    F.write('''
    model.component("comp1").geom("geom1").create("pt%i", "Point");
    model.component("comp1").geom("geom1").feature("pt%i").set("p", new String[]{"%i", "%i", "0"});    '''%(point_num,point_num,GPS[i,0],GPS[i,1]))
    point_num+=1
    run_count+=1
        
    if run_count==60:
        run_count=0
        run+=1
        F.write('''return model;
  }
  
  public static Model run%s(Model model) {
            '''%run)
F.write('''
    model.component("comp1").geom("geom1").create("rot1", "Rotate");
    model.component("comp1").geom("geom1").feature("rot1").set("pos", new String[]{"fault_x", "fault_y", "0"});
    model.component("comp1").geom("geom1").feature("rot1").set("rot", "strike");
    model.component("comp1").geom("geom1").feature("rot1").selection("input").set("wp1");
    
    model.component("comp1").mesh("mesh1").create("ftet1", "FreeTet");
    model.component("comp1").mesh("mesh1").create("ref1", "Refine");
    model.component("comp1").mesh("mesh1").create("ref2", "Refine");
    model.component("comp1").mesh("mesh1").feature("ftet1").selection().geom("geom1");
    model.component("comp1").mesh("mesh1").feature("ref2").selection().geom("geom1", 2);
    model.component("comp1").mesh("mesh1").feature("ref2").selection()
        .set(6''')
                                                                           
fault_index=7
for i in range(1,x_num*y_num):

    F.write(', %s'%fault_index)
    fault_index+=1

F.write(''');

    model.component("comp1").mesh("mesh1").feature("size").set("hauto", 3);
    model.component("comp1").mesh("mesh1").feature("size").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size").set("hmin", 150);
    model.component("comp1").mesh("mesh1").feature("ref1").set("boxcoord", true);
    model.component("comp1").mesh("mesh1").feature("ref1").set("xmax", "fault_y+fault_wid/2+1e3");
    model.component("comp1").mesh("mesh1").feature("ref1").set("xmin", "fault_y-fault_wid/2-1e3");
    model.component("comp1").mesh("mesh1").feature("ref1").set("ymax", "fault_y+fault_wid/2+1e3");
    model.component("comp1").mesh("mesh1").feature("ref1").set("ymin", "fault_y-fault_wid/2-1e3");
    model.component("comp1").mesh("mesh1").feature("ref1").set("zmin", "-fault_z-(fault_hei/2)-1e3");
    
    
    model.component("comp1").mesh("mesh1").feature("ref1").set("rmethod", "regular");
    model.component("comp1").mesh("mesh1").feature("ref1").set("numrefine", 2);
    model.component("comp1").mesh("mesh1").feature("ref2").set("rmethod", "regular");
    model.component("comp1").mesh("mesh1").feature("ref2").set("numrefine", 2);
    model.component("comp1").mesh("mesh1").run(); 
    model.component("comp1").geom("geom1").run("rot1");''')
    
    
def anomaly_maker(ano_size,ano_coor):
    F.write('''
    
    model.component("comp1").geom("geom1").create("blk2", "Block");
    model.component("comp1").geom("geom1").feature("blk2").set("size", new String[]{"%f", "1", "1"});
    model.component("comp1").geom("geom1").feature("blk2").setIndex("size", "%f", 1);
    model.component("comp1").geom("geom1").feature("blk2").setIndex("size", "%f", 2);
    model.component("comp1").geom("geom1").feature("blk2").set("base", "center");
    
    model.component("comp1").geom("geom1").feature("blk2").set("pos", new String[]{"%f", "0", "0"});
    model.component("comp1").geom("geom1").feature("blk2").setIndex("pos", "%f", 1);
    model.component("comp1").geom("geom1").feature("blk2").setIndex("pos", "%f", 2);
    model.component("comp1").geom("geom1").run("blk2");
    model.component("comp1").geom("geom1").runPre("fin");
        
    '''%(ano_size[0],ano_size[1],ano_size[2],ano_coor[0],ano_coor[1],ano_coor[2]))
    return
    
if anomaly == 'on':
    anomaly_maker(anomaly_size,anomaly_coor)
    
F.write('''
    
    model.label("green_mesh_Model.mph");

    model.component("comp1").geom("geom1").run();

    model.component("comp1").mesh("mesh1").run();
    model.component("comp1").mesh("mesh1").export().set("type", "nativeascii");
    model.component("comp1").mesh("mesh1").export("mesh.mphtxt");

    return model;
  }

  public static void main(String[] args) {
    Model model = run();''')
    
for i in range(2,run+1):
    F.write('''
    run%s(model);'''%i)
    
F.write('''
  }

}
''')
F.close()

############################################################################
##############################Run COMSOL####################################
os.system("comsol compile green_mesh.java")
os.system('comsol batch -inputfile green_mesh.class')
os.system('cp mesh.mphtxt pylith')
os.chdir('pylith')
############################################################################
############################################################################

def writing(name,x,y,z,dens,poisson,yong):
    
    G = yong/2/(1+poisson)
    K = yong/(3*(1-2*poisson))                       
    vp = sqrt((K+4/3*G)/dens)
    vs = sqrt(G/dens)
    
    matp = open('%s'%name,'w')
    matp.write('''#SPATIAL.ascii 1

SimpleDB {

  num-values = 3
  value-names =  density vs vp

  value-units =  kg/m**3  m/s  m/s

  num-locs = 1
 

  data-dim = 0

  space-dim = 3

  cs-data = cartesian {

    to-meters = 1.0

    space-dim = 3
  }
}
''')
    matp.write('%s %s %s %s %s %s \n\n' %(x,y,z,dens,vs,vp))
    matp.close()
    return()

writing('matprops.spatialdb',0,0,-domain_size[2]/2,domain_density,domain_nu,domain_E)

if anomaly=='on':
    writing('anomaly.spatialdb',anomaly_coor[0],anomaly_coor[1],anomaly_coor[2],anomaly_density,anomaly_nu,anomaly_E)


############################################################################
########################### mphtxt to mesh #################################

os.system('python convert.py %f %f %i %i %f %f %f %f %f %f %i %f %s'%(fault_wid, fault_hei, x_num, y_num, fault_coor[0], fault_coor[1], fault_coor[2], domain_size[0], domain_size[1], domain_size[2], dip_angle, -90-strike_angle, anomaly))

############################################################################
############################################################################


coor=np.loadtxt('./mesh_temp/coor.txt')
cell=np.loadtxt('./mesh_temp/cell.txt')
idd=np.loadtxt('./mesh_temp/id.txt')


F=open('../mesh.vtk','w')

F.write('''# vtk DataFile Version 2.0
Simplicial Mesh Example
ASCII
DATASET UNSTRUCTURED_GRID
POINTS %i double
'''%len(coor))

for i in range(len(coor)):
    F.write('%f %f %f \n'%(coor[i,1],coor[i,2],coor[i,3]))
    
F.write('''
CELLS %i %i
'''%(len(cell),5*len(cell)))

for i in range(len(cell)):
    F.write('4  %i %i %i %i \n'%(cell[i,1],cell[i,2],cell[i,3],cell[i,4]))
    
F.write('''
CELL_TYPES %i
'''%(len(cell)))

for i in range(len(cell)):
    F.write('10 \n')
    
F.write('''CELL_DATA %i
SCALARS density double 1
LOOKUP_TABLE default
'''%len(cell))
for i in range(len(cell)):
    if idd[i,1] == 1:
        F.write('%f\n'%domain_density)
    elif idd[i,1] == 2:
        F.write('%f\n'%anomaly_density)

F.write('''SCALARS nu double 1
LOOKUP_TABLE default
''')

for i in range(len(cell)):
    if idd[i,1] == 1:
        F.write('%f\n'%domain_nu)
    elif idd[i,1] == 2:
        F.write('%f\n'%anomaly_nu)

F.write('''SCALARS E double 1
LOOKUP_TABLE default
''')

for i in range(len(cell)):
    if idd[i,1] == 1:
        F.write('%f\n'%domain_E)
    elif idd[i,1] == 2:
        F.write('%f\n'%anomaly_E)
    else:
        print('error error error error error')
F.close()
    
########################## Run PyLith ######################################
    
for i in range(x_num*y_num):
    os.system("pylith fault_%i.cfg --nodes=8"%(i+1))

##################### Make Green's function ################################
    
def make_Green(GREENFN):
    tong=np.array(GREENFN,dtype=np.float64)
    tong=tong.reshape(x_num*y_num,-1)
    greenfn=tong.T
    return(greenfn)

os.chdir("output")

tol = 1

GPS_X = []
GPS_Y = []
GPS_Z = []

for i in range(1,x_num*y_num+1): ####  
    
    F=open("dislocation%i_t1000000.vtk"%(i),"r")
    lines = F.readlines()
    F.close()
    
    for j in range(len(lines)):
        if "POINTS" in lines[j]:
            PS = j+1
        elif "CELLS" in lines[j]:
            PS_END = j
        elif "VECTORS" in lines[j]:
            VS = j+1
            
    POINT = lines[PS:PS_END]
    POINT = np.array([l.split(' ') for l in POINT],dtype=np.float64)
    
    DISPL = lines[VS:]
    DISPL = np.array([l.split(' ') for l in DISPL],dtype=np.float64)
    
    
    surface=[]
    displace=[]
    for zero in range(len(POINT)):            ####Find Surface Coordinate
        if POINT[zero,2] > -1e-3:
            surface.append(POINT[zero,0:3])
            displace.append(DISPL[zero,0:3])
    POINT=np.array(surface)
    DISPL=np.array(displace)

    observe_point=0
    for k in range(len(GPS)):
        for j in range(len(POINT)):
            x=POINT[j,0]
            y=POINT[j,1]
            if GPS[k,0]-tol <x< GPS[k,0]+tol and GPS[k,1]-tol<y<GPS[k,1]+tol :
                GPS_X.append(DISPL[j,0])
                GPS_Y.append(DISPL[j,1])
                GPS_Z.append(DISPL[j,2])
                observe_point+=1
                break
        print(observe_point)
             
    
grfn_X = make_Green(GPS_X)
grfn_Y = make_Green(GPS_Y)
grfn_Z = make_Green(GPS_Z)

grfn_hap = np.vstack((np.vstack((grfn_X,grfn_Y)),grfn_Z))

os.chdir(home)
np.savetxt('green_function.dat',grfn_hap,delimiter = ' ')




    
    
    
    
    
    
