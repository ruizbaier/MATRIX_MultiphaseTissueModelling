lc=0.25;
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {0, 1, 0, lc};


Line(11) = {1, 2};
Circle(12) = {2, 1, 3};
Line(13) = {3, 1};

Line Loop(14) = {11, 12, 13};
Plane Surface(15) = {14};
Physical Surface(20) = {15};

//bottom 
Physical Line(31) = {11};
//circle
Physical Line(32) = {12};
//left 
Physical Line(33) = {13};

Mesh.CharacteristicLengthFactor = 0.04;
Mesh.ScalingFactor = 1.0;
Mesh.Smoothing = 10;
Mesh.Algorithm = 2;