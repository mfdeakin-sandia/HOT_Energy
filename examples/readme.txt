HOT_Energy environment variable should be defined to the local repo clone.
e.g. /Users/samitch/Documents/repos/PrimalDual/HOT_Energy
If you run from within the debugger, you'll need to set
Product > Scheme > Edit Scheme > Run in the sidebar > Arguments
-- 

Example meshes have three files

*_points.txt are the x,y coordinates of each of the points, numbered 1..n
*_weights.txt are the weights of each of the corresponding points 1..n
*_triangles.txt are the indices of the points in each triangle 1..t, ccw order

e.g.

0 0
1 0
0 1

0
0
0

1 2 3