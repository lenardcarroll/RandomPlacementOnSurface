# RandomPlacementOnSurface

This is a script I originally wrote to just place random oxygen atoms on top of a Au(221) surface. I have reworked it and it might fit your needs. How it works:

Any molecule you want added to the surface, place it into a file with extension .txt. This file must have a format of:
```
N
A x1 y1 z1
B x2 y2 z2
C x3 y3 z3
...
Z xk yk zk
```

where ```N``` is the number of copies of the molecule in this file that should be placed on the surface, ```A-Z``` is the atoms that make up the molecule and next to it, is its coordinates. Make sure that the atom that will be mostly interacting with the surface is zeroed, i.e., it is at the origin, with no other atom having a lower z-value.

When using the script, it is a general form of:

```
python -i RandomPlacementOnSurface.py -in <STRUCTURE_FILE_WITH_SURFACE> -out <NAME_OF_OUTPUT_FILE> -plot <Y or N> -mdist <REAL_VALUE_MIN_DISTANCE_BETWEEN_SURFACE_ATOMS_AND_ADSORBED_ATOMS> -xdist <REAL_VALUE_MAX_DISTANCE_BETWEEN_SURFACE_ATOM_AND_ADSORBED_ATOM> -adist <REAL_VALUE_MIN_DISTANCE_BETWEEN_ADSORBED_ATOMS>
```

with ```-plot``` you can decide if you want to make a simple plot of the surface (Y) or not (N). ```-mdist``` is purely the minimum distance between the adsorbed atom and the surface atoms, while ```-xdist``` is the maximum distance, and ```-adist``` is the minimum distance between adsorbed atoms (these being the central atoms). 

The script works my reading in the structure file, but sure it only contains the surface of the structure. Next, it finds the outermost points of this surface via a convex hull algorithmn. These points are used to create a polygon, with this polygon split into a bunch of triangles. Points are then placed randomly in a random triangle using the equation:

<x4,y4,z4> = a1<x1,y1,z1> + a2<x2,y2,z2>

where <x3,y3,z3> = <0,0,0>, and a1,a2 ∈ [0,1], and a1+a2<=1.

On the random points generated, the moelcules are placed and the final output structure generated and saved. Below is a simple plot and a 3D image.

![Polygon Surface with Adsorbed Points](https://raw.githubusercontent.com/lenardcarroll/RandomPlacementOnSurface/main/demonstration3.jpg "Polygon Surface with Triangles")
![Polygon Surface with Adsorbed Points](https://raw.githubusercontent.com/lenardcarroll/RandomPlacementOnSurface/main/demonstration4.jpg "Polygon Surface with Triangles and Adsorbed Points")
![Polygon Surface with Adsorbed Points](https://raw.githubusercontent.com/lenardcarroll/RandomPlacementOnSurface/main/demonstration2.jpg "Polygon Surface with Adsorbed Points")
![Surface with adsorbed molecules](https://raw.githubusercontent.com/lenardcarroll/RandomPlacementOnSurface/main/Example13D.png "Surface with adsorbed molecules")
![Surface with adsorbed molecules - Part2](https://raw.githubusercontent.com/lenardcarroll/RandomPlacementOnSurface/main/Example23D.png "Surface with adsorbed molecules - Different angle")
