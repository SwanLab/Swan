function disp = CohesiveSeparationComputer(u,cohesiveMesh)
% u es una LagrangianFunction definida sobre una cohesive mesh

%Necessito crear una malla 1D sobre els punts mitjos dels 

coord = cohesiveMesh.mesh.coord;
nMidNodes = size(cohesiveMesh.pairsMatrix,1);

midCoord= (coord(cohesiveMesh.listNodeCohesive,:)+ ...
    coord(cohesiveMesh.pairsMatrix(:,2),:))/2;

    midCoord = midCoord(:,1);



midConnec =  [(1:nMidNodes-1)' (2:nMidNodes)'];

s.connec = midConnec;
s.coord = midCoord;

        

midMesh = Mesh.create(s);

sep = LagrangianFunction.create(midMesh,2,'P1');

% posar setFValues al disp
% he de mirar com es defineixen els separations


disp =1;

end