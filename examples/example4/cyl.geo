lc = 0.02;
r =0.915;
t = 0.21;
nl = 9;

Point(1) = {0.,  0., -t/2., lc};
Point(2) = {r,   0., -t/2., lc};
Point(3) = {-r,  0., -t/2., lc};
Point(4) = {0.,  r,  -t/2., lc};
Point(5) = {0., -r,  -t/2., lc};
Circle(1) = {2,1,4};
Circle(2) = {4,1,3};
Circle(3) = {3,1,5};
Circle(4) = {5,1,2};

Line Loop (1) = {1,2,3,4};

Plane Surface(1) = {1};

Extrude {0,0,t} {Surface{1}; Layers{nl};}
Physical Volume(1)={1};
