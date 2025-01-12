// Gmsh project created on Wed Jan 08 17:01:54 2025
Mesh.ElementOder =1;
Mesh.Algorithm = 5;

mesh_edge=40;
//+
Point(1) = {0, 0, 0};
//+
Point(2) = {1, 0, 0};
//+
Point(3) = {1, 1, 0};
//+
Point(4) = {0, 1, 0};
//+

//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+

Transfinite Line{1, 2, 3, 4} = mesh_edge+1 Using Progression 1;
Transfinite Surface{1};
Recombine Surface{1};//+
Physical Curve("left", 5) = {1};
//+
Physical Curve("top", 6) = {4};
//+
Physical Curve("right", 7) = {3};
//+
Physical Curve("bottom", 8) = {2};
//+
Physical Surface("surface", 9) = {1};
