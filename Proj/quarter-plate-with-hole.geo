// Gmsh project created on Wed Jan 08 21:15:36 2025
R = 0.5;
L = 2.0;

Point(1) = {L, 0, 0};
Point(2) = {L, L, 0};
Point(3) = {0, L, 0};
Point(4) = {0, 0, 0};
Point(5) = {R, 0, 0};
Point(6) = {0, R, 0};
Point(7) = { Cos(Pi/4) * R,  Sin(Pi/4) * R, 0};

Circle(1) = {5, 4, 7};
Circle(2) = {7, 4, 6};

Line(3) = {6, 3};
Line(4) = {3, 2};
Line(5) = {2, 1};
Line(6) = {1, 5};
Line(7) = {2, 7};

Curve Loop(1) = {4, 7, 2, 3};
Plane Surface(1) = {1};

Curve Loop(2) = {7, -1, -6, -5};
Plane Surface(2) = {2};

Transfinite Line{1, 2, 3, 4, 5, 6, 7} = 20 +1;

Transfinite Surface{1};
Transfinite Surface{2};

Recombine Surface{1};
Recombine Surface{2};

Mesh.ElementOrder = 1;
Mesh.Algorithm = 8;//+

Physical Curve("circle", 8) = {2, 1};
//+
Physical Curve("bottome", 9) = {6};
//+
Physical Curve("right", 10) = {5};
//+
Physical Curve("top", 11) = {4};
//+
Physical Curve("left", 12) = {3};
//+
Physical Curve("Oblique", 13) = {7};
//+
Physical Surface("surface", 14) = {1, 2};
