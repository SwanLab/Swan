function [L,D] = aero_forces(nodespresscyl,pressureFun,mesh)

nodesCyl    = nodespresscyl; 
presCylVals = pressureFun.fValues(nodesCyl,1);
%xCyl        = mesh.coord(nodesCyl,1);
%yCyl        = mesh.coord(nodesCyl,2);
mesh.computeEdges();
e  = mesh.edges.nodesInEdges; %Get the List of All Edges
bE = ismember(e,nodesCyl); 
bE = find(prod(bE,2));   % Find Edges Where Both Nodes Belong to the Cylinder
connec = e(bE,:); %stores the edges that form the cylinder's boundary.
ss.coord    = mesh.coord;
ss.connec   = connec;
ss.kFace    = -1;
bMesh       = Mesh.create(ss);
bMesh       = bMesh.computeCanonicalMesh();
presCyl     = LagrangianFunction.create(bMesh,1,pressureFun.order); 
presCyl.fValues = presCylVals;

presCyl.plot()

normal_vectors = zeros(bMesh.nelem,bMesh.ndim); % define les matrius de vector normals i la longitud d'element.
length_element = zeros(bMesh.nelem,1);

centroid = mean(bMesh.coord); % calcula el centroide
central_points = (bMesh.coord(bMesh.connec(:,1),:)+bMesh.coord(bMesh.connec(:,2),:))/2; % punt central de cada edge
%ref_vect = central_points - centroid;

cont =1;

for iE = 1:bMesh.nelem
    node1 = bMesh.coord(bMesh.connec(iE,1),:); %  two endpoints of the edge.
    node2 = bMesh.coord(bMesh.connec(iE,2),:);

    %if node1(1)<= 5
    nvect = (node2-node1)/(abs(norm(node2-node1))); % busca el vector tangent al edge
    nvect = -nvect * [0 -1;1 0]; % el gira -90º per obtenir el vector normal
%     if dot(ref_vect(iE,:),nvect)<0 %No cal
%         nvect = -nvect;
%     end
    normal_vectors(cont,:) = nvect;
    length_element(cont) = abs(norm(node1-node2));

    cont = cont +1;
    %end

end

nx = LagrangianFunction.create(bMesh,1,'P0');%Vectors normals
ny = LagrangianFunction.create(bMesh,1,'P0');
nx.fValues = normal_vectors(:,1);
ny.fValues = normal_vectors(:,2);
sss.operation = @(x) -presCyl.evaluate(x).*nx.evaluate(x);
pNx           = DomainFunction(sss);
D             = Integrator.compute(pNx,bMesh,2);
sss.operation = @(x) -presCyl.evaluate(x).*ny.evaluate(x);
pNy           = DomainFunction(sss);
L             = Integrator.compute(pNy,bMesh,2);

quiver(central_points(:,1),central_points(:,2),normal_vectors(:,1),normal_vectors(:,2)) %Plot the vectors
hold on
quiver(centroid(1,1),centroid(1,2),D,0);
hold on
quiver(centroid(1,1),centroid(1,2),0,L);
hold on
bMesh.plot() %Plot mesh points



end



































