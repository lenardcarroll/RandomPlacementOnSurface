#All modules/packages needed
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import random
from pathlib import Path
import math
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

#Arguments Cheat/Help Sheet
parser = argparse.ArgumentParser()
parser.add_argument("-in", "--input", dest = "input", default = "input.xyz", help="Name of the input structure file, which contains the coordinates.")
parser.add_argument("-mdist", "--mindist", dest = "mindist", default = "2.15", help="Specify the minimum allowed distance between surface atoms and adsorbed atoms")
parser.add_argument("-xdist", "--maxdist", dest = "maxdist", default = "2.50", help="Specify the maximum allowed distance between surface atoms and adsorbed atoms")
parser.add_argument("-adist", "--adist", dest = "adist", default = "4.25", help="Specify the minimum allowed distance between adsorbed atoms")
parser.add_argument("-out", "--output", dest = "output", default = "output.xyz", help="Choose the file you want the final coordinates to be printed out to.")
parser.add_argument("-plot", "--plot", dest = "plot", default = "N", help="If Y is chosen, a simple plot is made of the surface and the random points. Otherwise (N) none is made.")

#parser.add_argument("-conf", "--confim", dest = "confirm", default = "N", help="If Y, the script will confirm that no points lie outside the convex, if it does, it will let you know. N is used to not do the check.")
args = parser.parse_args()

how_many_mol = []
num_lines = []

'''
The script reads in every .txt file in a directory as it assumes it contains information on the molecules you want to place on the surface. Make sure to place
only the molecule .txt files in the directory containing the python script. Furthermore, the .txt file must be of a specific form, with that form looking like the following:

A
C x1 y1 z1
O x2 y2 z2
...
C xn yn zn

Where A is the number of copies of the molecule in the current .txt file that should be placed on the surface, and the rest are coordinates of the molecule. Make sure that
the atom you want to mainly interact with the surface has been set at the origin, particularly, it is the lowest coordinate in the molecule. 
'''
for p in Path('.').glob('*.txt'):
    with open(p.name) as f:
        lines = f.read()
        first = lines.split('\n', 1)[0]
    how_many_mol.append(int(first))
    num_lines.append(len(list(open(p.name))))
    if len(how_many_mol) == 1:
        df = pd.read_csv(p.name, skiprows=1, names=['Atom', 'X', 'Y', 'Z'], sep="\s+", engine='python')
    else:
        df2 = pd.read_csv(p.name, skiprows=1, names=['Atom', 'X', 'Y', 'Z'], sep="\s+", engine='python')
        df = pd.concat([df, df2], axis=0)

#Split all the loaded in molecules (to be added onto the surface) into separate parts.
sum_lines = 0
range_info = []
for i in range(len(num_lines)):
    range_info.append([sum_lines,sum_lines+num_lines[i]-2])
    sum_lines += num_lines[i] - 1
frames = [ df.iloc[i[0]:i[1]+1].copy() for i in range_info ]

#Determines how many times a specific molecule should be placed on the surface
How_Many_Systems = []
for i in range(len(how_many_mol)):
    for j in range(how_many_mol[i]):
        How_Many_Systems.append(i)

#Load the surface in
df = pd.read_csv(args.input, skiprows=2, names=['Atom', 'X', 'Y', 'Z'], sep="\s+", engine='python')
df = df.sort_values(by ='X', ascending = True)

coordinates = []
for i in df.index:
    coordinates.append([str(df['Atom'].iloc[i]), df['X'].iloc[i], df['Y'].iloc[i], df['Z'].iloc[i]])
    
#Added sorted X and Y coordinates to list points
points = []
for i in range(len(df)):
    points.append([df['X'].iloc[i],df['Y'].iloc[i]])

#Add first point into list convex_points
convex_points = [[points[0][0],points[0][1]]]

#Give the first point an angle 90
angles = [90]

#Algorithm starts, check https://github.com/lenardcarroll/myConvexHull.py/blob/main/README.md for an explanation
j = 0
while j>-1:
    altered_angles = []
    k_points = []
    unaltered_angles = []
    i = convex_points[len(convex_points)-1]
    for k in range(len(points)):
        if i!=points[k]:
            altx1y1 = [i[0] - i[0],i[1] - i[1]]
            altx2y2 = [points[k][0] - i[0],points[k][1] - i[1]]
            altx2y2 = altx2y2/np.linalg.norm(altx2y2)
            ang = np.arctan2(altx2y2[1],altx2y2[0]) - np.arctan2(altx1y1[1],altx1y1[0])
            altang = ang - angles[j] + 2*np.pi
            if altang > 2*np.pi:
                altang = altang - 2*np.pi
            altered_angles.append(altang)
            unaltered_angles.append(ang)
            k_points.append(k)
    max_altered_angle = max(altered_angles)
    convex_points_max = points[k_points[altered_angles.index(max_altered_angle)]]
    convex_points.append(convex_points_max)
    convex_points_unaltered_angle = unaltered_angles[altered_angles.index(max_altered_angle)]
    angles.append(convex_points_unaltered_angle)

    if convex_points_max == convex_points[0]:
        j=-1
    if j > 0:
        points.remove(convex_points_max)
        j+=1
    elif j == 0:
        j+=1

#Construct a polygon out of the convex points
polygonsurface = []
for i in convex_points:
    polygonsurface.append((i[0],i[1]))
polygon = Polygon(polygonsurface)

#Find the midpoint of the polygon
sumMidPointx = 0
sumMidPointy = 0
for i in convex_points:
    sumMidPointx += i[0]
    sumMidPointy += i[1]   
MidPoint = [sumMidPointx/len(convex_points),sumMidPointy/len(convex_points)]

#Split the polygon into a bunch of triangles using the polygon vertices
Triangles = []
for i in range(1,len(convex_points)):
    Triangles.append([convex_points[i-1],convex_points[i],MidPoint])

i=0

OO_distance = float(args.adist)
O_atoms = []

#function to calculate the distance between multiple vectors
def dist2(a,b,c):
    if np.abs(np.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)-np.sqrt((a[0] - c[0])**2 + (a[1] - c[1])**2)-np.sqrt((c[0] - b[0])**2 + (c[1] - b[1])**2))<0.001:
        return True
'''
This process/algorithm works as follows:
Randomly select a triangle from the polygon and then take this move this triangle such that a point on the triangle (the midpoint from the polygon) is at the origin. 
Then use the equation: <x4,y4,z4> = a1*<x1,y1,z1> + a2<x2,y2,z2>
to find a point inside the triangle, where a1 and a2 are real values between 0 and 1, and a1 + a2 = 1. This point is then shifted back to where the triangle used to be in
the polygon. Next, the algorithm makes sure the point is either in the polygon or on the lines connecting the vertices. 
'''
while i<len(How_Many_Systems):
    TriangleChoice = random.choice(range(len(Triangles)))
    j = Triangles[TriangleChoice]
    x1y1 = [j[0][0]-j[2][0],j[0][1]-j[2][1]]
    x2y2 = [j[1][0]-j[2][0],j[1][1]-j[2][1]]
    MidPoint = [0,0]
    r = 0
    while r>-1:
        a1 = random.uniform(0,1)
        a2 = random.uniform(0,1)
        if (a1+a2)>1:
            r+=1
        else:
            r=-1
    newxy = [a1*x1y1[0]+a2*x2y2[0],a1*x1y1[1]+a2*x2y2[1]]
    newxy = [newxy[0]+j[2][0],newxy[1]+j[2][1]]
    if polygon.contains(Point(newxy[0],newxy[1])) == True:
        z_preliminary = random.uniform(2.15, 2.50)
        atom_distances = []
        Error_list = []
        for k in range(len(coordinates)):
            atom_distances.append(np.sqrt((coordinates[k][1]-newxy[0])**2+(coordinates[k][2]-newxy[1])**2))
        closest_member = min(atom_distances)
        if closest_member > z_preliminary:
            i+=0
            continue
        atom_number = 0
        for k in range(len(coordinates)):
            if closest_member == atom_distances[k]:
                atom_number = k
        z_value = np.sqrt((z_preliminary)**2 - (coordinates[atom_number][1]-newxy[0])**2 - (coordinates[atom_number][2]-newxy[1])**2) + coordinates[atom_number][3]
        if math.isnan(z_value) == True:
            i+=0
            continue
        for k in range(len(coordinates)):
            value = np.sqrt((coordinates[k][1]-newxy[0])**2+(coordinates[k][2]-newxy[1])**2+(coordinates[k][3]-z_value)**2)
            if value<2.15:
                Error_list.append("Error")
        if len(Error_list)>0:
            i+=0
            continue
        if len(O_atoms)>0:
            for k in range(len(O_atoms)):
                value = (np.sqrt((O_atoms[k][1]-newxy[0])**2+(O_atoms[k][2]-newxy[1])**2+(O_atoms[k][3]-z_value)**2))
                if value<OO_distance:
                    Error_list.append("Error")
        if len(Error_list)>0:
            i+=0
            continue
        else:
            O_atoms.append(['O',newxy[0],newxy[1],z_value])
            i += 1
    else:
        convex_close = []
        for j in convex_points:
            dist = np.abs((newxy[0]-j[0])**2+(newxy[1]-j[1])**2)
            if dist<0.001:
                convex_close.append(dist)
        if len(convex_close) == 0:
            dotp = []
            for k in range(len(convex_points)):
                for l in range(len(convex_points)):
                    if k!=l:
                        if dist2(convex_points[k],convex_points[l],[newxy[0],newxy[1]]) == True:
                            dotp.append(1)
            if len(dotp)==0:
                i+=0
                continue
            else:
                z_preliminary = random.uniform(float(args.mindist), float(args.maxdist))
                atom_distances = []
                Error_list = []
                for k in range(len(coordinates)):
                    atom_distances.append(np.sqrt((coordinates[k][1]-newxy[0])**2+(coordinates[k][2]-newxy[1])**2))
                closest_member = min(atom_distances)
                if closest_member > z_preliminary:
                    i+=0
                    continue
                atom_number = 0
                for k in range(len(coordinates)):
                    if closest_member == atom_distances[k]:
                        atom_number = k
                z_value = np.sqrt((z_preliminary)**2 - (coordinates[atom_number][1]-newxy[0])**2 - (coordinates[atom_number][2]-newxy[1])**2) + coordinates[atom_number][3]
                if math.isnan(z_value) == True:
                    i+=0
                    continue
                for k in range(len(coordinates)):
                    value = np.sqrt((coordinates[k][1]-newxy[0])**2+(coordinates[k][2]-newxy[1])**2+(coordinates[k][3]-z_value)**2)
                    if value<float(args.mindist):
                        Error_list.append("Error")
                if len(Error_list)>0:
                    i+=0
                    continue
                if len(O_atoms)>0:
                    for k in range(len(O_atoms)):
                        value = (np.sqrt((O_atoms[k][1]-newxy[0])**2+(O_atoms[k][2]-newxy[1])**2+(O_atoms[k][3]-z_value)**2))
                        if value<OO_distance:
                            Error_list.append("Error")
                if len(Error_list)>0:
                    i+=0
                    continue
                else:
                    O_atoms.append(['O',newxy[0],newxy[1],z_value])
                    i += 1
        else:
            z_preliminary = random.uniform(float(args.mindist), float(args.maxdist))
            atom_distances = []
            Error_list = []
            for k in range(len(coordinates)):
                atom_distances.append(np.sqrt((coordinates[k][1]-newxy[0])**2+(coordinates[k][2]-newxy[1])**2))
            closest_member = min(atom_distances)
            if closest_member > z_preliminary:
                i+=0
                continue
            atom_number = 0
            for k in range(len(coordinates)):
                if closest_member == atom_distances[k]:
                    atom_number = k
            z_value = np.sqrt((z_preliminary)**2 - (coordinates[atom_number][1]-newxy[0])**2 - (coordinates[atom_number][2]-newxy[1])**2) + coordinates[atom_number][3]
            if math.isnan(z_value) == True:
                i+=0
                continue
            for k in range(len(coordinates)):
                value = np.sqrt((coordinates[k][1]-newxy[0])**2+(coordinates[k][2]-newxy[1])**2+(coordinates[k][3]-z_value)**2)
                if value<float(args.mindist):
                    Error_list.append("Error")
            if len(Error_list)>0:
                i+=0
                continue
            if len(O_atoms)>0:
                for k in range(len(O_atoms)):
                    value = (np.sqrt((O_atoms[k][1]-newxy[0])**2+(O_atoms[k][2]-newxy[1])**2+(O_atoms[k][3]-z_value)**2))
                    if value<OO_distance:
                        Error_list.append("Error")
            if len(Error_list)>0:
                i+=0
                continue
            else:
                O_atoms.append(['O',newxy[0],newxy[1],z_value])
                i += 1

new_O_atoms = []

'''
Place the molecules from the .txt files in the positions randomly generated on the surface. 
'''

total_mol = []
for i in range(len(How_Many_Systems)):
    total_mol.append(i)
for i in range(len(O_atoms)):
    x = random.choice(total_mol)
    for j in frames[How_Many_Systems[x]].index:
        X = frames[How_Many_Systems[x]]['X'].iloc[j] + O_atoms[i][1]
        Y = frames[How_Many_Systems[x]]['Y'].iloc[j]  + O_atoms[i][2]
        Z = frames[How_Many_Systems[x]]['Z'].iloc[j]  + O_atoms[i][3]
        new_O_atoms.append([frames[How_Many_Systems[x]]['Atom'].iloc[j],X,Y,Z])
    total_mol.remove(x)

f = open(args.output, 'w')
print(len(coordinates)+len(new_O_atoms),file=f)
print("Surface Generated",file=f)
for k in coordinates:
    print(k[0],k[1],k[2],k[3],file=f)
for i in new_O_atoms:
    print(i[0], i[1], i[2], i[3],file=f)
f.close()

if args.plot == 'Y':
    X1 = []
    Y1= []
    for i in convex_points:
        X1.append(i[0])
        Y1.append(i[1])
    plt.scatter(X1,Y1,color='blue',lw=3)
    plt.plot(X1,Y1,color='blue',lw=3)
    for i in O_atoms:
        plt.scatter(i[1],i[2],color='red',edgecolors='black',s=250,lw=2)
    plt.gcf().set_size_inches(19.2, 14.4)
    plt.grid(False)
    plt.xlabel("x-coordinates",fontsize=20)
    plt.ylabel("y-coordinates",fontsize=20)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.show()
