function createRegularHexagonMesh()
meshfilename = 'test2d_micro_joel_PRUEBA_hexagon.m';
Data_prb = {'''TRIANGLE''','''SI''','''2D''','''Plane_Stress''','''ELASTIC''','''MICRO'''};
ndim = 2;
c = 3;
div = 2; % Divisiones por lado del hexagono
% Tres divisiones equivalen a cuatro nodos por lado 
ydiv = div;
xdiv = div;
nnodes = 3*(div+1)*div+1;
coord = zeros(nnodes,ndim); %Se debe cambiar el número de prueba
row = 1;
for jdiv = 0:ydiv
    for idiv = 0:xdiv
        coord(row,1) = coord(row,1) + 0.5*c + c/div*idiv - 0.5*c/div*jdiv;
        coord(row,2) = coord(row,2) + 2*c*0.866 - 0.866*c/div*jdiv;
        row = row+1;
    end
    xdiv = xdiv+1;
end
xdiv = xdiv-1;
for jdiv = 1:ydiv
    xdiv = xdiv-1;
    for idiv = 0:xdiv
        coord(row,1) = coord(row,1) + c/div*idiv + 0.5*c/div*jdiv;
        coord(row,2) = coord(row,2) + c*0.866 - 0.866*c/div*jdiv;
        row = row+1;
    end
end
connec = delaunay(coord);
s.coord = coord;
s.connec = connec;
m = Mesh(s);
m.plot();
master_slave = computeMasterSlaveNodes(coord,div,c);
end

