% Generate mesh
    coord = [0,0; %1
             1,0; %2
             2,0; %3
             0,1; %4
             1,1; %5
             2,1; %6
             0,2; %7
             1,2; %8
             2,2;]; %9
    
    connec = [1 2 5 4;
              2 3 6 5;
              4 5 8 7;
              5 6 9 8];
    
    s.coord  = coord;
    s.connec = connec;
    m = Mesh(s);   

% Quad mesh to tri
mtri = m.convertToTriangleMesh();