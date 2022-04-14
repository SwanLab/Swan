function SquaredMeshCreator()
    [dim,divUnit,c,theta] = obtainInitialData();
    nsides = obtainPolygonSides(c,theta);
    [coord,vertCoord,boundary,boundNodes,div] = initializeVariables(dim,divUnit,nsides,c); 
    vertCoord = computeVertCoord(vertCoord,c,theta,nsides);
    boundary = computeBoundaryCoord(boundary,vertCoord,c,theta,nsides,div); 
    coord = computeMeshCoord(nsides,vertCoord,divUnit,c,boundary,boundNodes,coord,div);
    connec = computeConnectivities(coord);
    plotCoordinates(coord,connec);
    masterSlaveIndex = obtainMasterSlaveNodes(vertCoord,boundary,nsides,div,dim,boundNodes); 
    vertIndex(:,1) = 1:nsides;
    plotVertices(vertIndex,coord);
    plotMasterSlaveNodes(masterSlaveIndex,coord);
    writeFEMreadingfunction(coord,connec,masterSlaveIndex,'Hexagon5x5x5.m', vertCoord);
end

function  [dim,divUnit,c,theta] = obtainInitialData()
% Datos de entrada del programa. COMPLETAMENTE GENERAL
    dim = 2;
    divUnit = 5; %Divisions/length of the side
    c = [1,1,1];
    theta = [0,60,120];
end

function nsides = obtainPolygonSides(c,theta)
% Obtención de nº de lados y filtro. COMPLETAMENTE GENERAL
    if length(c) == length(theta)
        nsides = 2*length(c);
    else
        cprintf('red','CRYTICAL ERROR. Vectors c and theta must have the same length\n');
    end
end

function [coord,vertCoord,boundary,boundNodes,div] = initializeVariables(dim,divUnit,nsides,c)
% Inicialización de variables. GENERAL a los casos deseados
div = divUnit*c;
    switch nsides
        case 4
            divA = div(1);
            divB = div(2);
            boundNodes = nsides*(1+1/2*(divA+divB-2));
            nnodes = boundNodes+(divA-1)*(divB-1);
        case 6
            divA = div(1);
            divB = div(2);
            divC = div(3);
            boundNodes = nsides*(1+1/3*(divA+divB+divC-3));
            nnodes = 0;
            sideNodes = div-1;
            while max(sideNodes) >= 1
                for i = 1:length(sideNodes)
                    if sideNodes(i) <= 0
                        sideNodes(i) = 0;
                    end
                end
                nnodes = nnodes+nsides+sum(2*sideNodes);
                sideNodes = sideNodes-1;
            end
            nnodes = nnodes+nsides+1;
    end
    boundNodes = round(boundNodes);
    vertCoord = zeros(nsides,dim);
    boundary = zeros(boundNodes,dim);
    coord = [];
end

function vertCoord = computeVertCoord(vertCoord,c,theta,nsides)
% Obtención de los vertices. COMPLETAMENTE GENERAL
    c0 = [0,0];
    for iMaster = 1:nsides/2
        pos = computeThePosition(c0,c(iMaster),theta(iMaster));
        vertCoord(iMaster+1,:) = vertCoord(iMaster+1,:)+pos;
        c0 = pos;
    end
    for iSlave = 1:nsides/2
        pos = computeThePosition(c0,c(iSlave),theta(iSlave)+180);
        if iSlave == nsides/2
            if vertCoord(1,:) ~= pos
                cprintf('red','CRYTICAL ERROR. Vertices computed wrongly \n');
            end
        else
            vertCoord(iMaster+iSlave+1,:) = vertCoord(iMaster+iSlave+1,:)+pos;
            c0 = pos;
        end
    end
end

function boundary = computeBoundaryCoord(boundary,vertCoord,c,theta,nsides,div)
% Obtención de las coord de la boundary. COMPLETAMENTE GENERAL
    boundary(1:nsides,:) = boundary(1:nsides,:)+vertCoord;
    cont = nsides+1;
    for iMaster = 1:nsides/2
        c0 = vertCoord(iMaster,:);
        for iDiv = 1:div(iMaster)-1
            pos = computeThePosition(c0,c(iMaster)/div(iMaster),theta(iMaster));
            boundary(cont,:) = boundary(cont,:)+pos;
            cont = cont+1;
            c0 = pos;
        end
    end
    for iSlave = 1:nsides/2
        c0 = vertCoord(iMaster+iSlave,:);
        for iDiv = 1:div(iSlave)-1
            pos = computeThePosition(c0,c(iSlave)/div(iSlave),theta(iSlave)+180);
            boundary(cont,:) = boundary(cont,:)+pos;
            cont = cont+1;
            c0 = pos;
        end
    end
end

function pos = computeThePosition(c0,c,theta)
    pos = c0+c.*[cosd(theta) sind(theta)];
end

function coord = computeMeshCoord(nsides,vertCoord,divUnit,c,boundary,boundNodes,coord,div)
% Obtención de las coordenadas dentro de la boundary. GENERAL para los
% casos considerados
coord(1:boundNodes,:) = boundary;
intNode = boundNodes+1;
    switch nsides
        case 4
        % Compute coords by intersections
            vA = vertCoord(2,:)-vertCoord(1,:);
            mA = norm(vA);
            vA = vA/mA;
            vB = vertCoord(3,:)-vertCoord(2,:);
            mB = norm(vB);
            vB = vB/mB;
            nodesX = divUnit*c(1)-1;
            nodesY = divUnit*c(2)-1;
            for jNodes = 1:nodesY
                pB = boundary(boundNodes+1-jNodes,:);
                for iNodes = 1:nodesX
                    pA = boundary(nsides+iNodes,:);
                    if vA(1) == 0
                        x = pB(1);
                        y = (x-pA(1))*vB(2)/vB(1)+pA(2);
                    elseif vB(1) == 0
                        x = pA(1);
                        y = (x-pB(1))*vA(2)/vA(1)+pB(2);
                    else
                        x = (pB(2)-pA(2)+pA(1)*vB(2)/vB(1)+pB(1)*vA(2)/vA(1))/(vB(2)/vB(1)-vA(2)/vA(1));
                        y = (x-pB(1))*vA(2)/vA(1)+pB(2);
                    end
                    coord(intNode,:) = [x y];
                    intNode = intNode+1;
                end
            end

        case 6
        % Compute coords by diagonals
        % Sitúa el nodo central
        vA = vertCoord(4,:)-vertCoord(1,:);
        pA = vertCoord(1,:);
        centerVec = vA/2;
        O = pA+centerVec;
        coord(intNode,:) = O;
        % Cálculo de las divisiones a aplicar por recta
        diagDiv = max(div);
        div = div-1;
        intNode = intNode+1;
        % Aplicación de las divisiones por cada semirecta
        for iDiv = 1:diagDiv-1
            for iDiag = 1:nsides
                diagA = O-vertCoord(iDiag,:);
                vecDiv = iDiv*diagA/diagDiv;
                pos = vertCoord(iDiag,:)+vecDiv;
                coord(intNode,:) = pos;
                intNode = intNode+1;
            end
            % Aplicar aquí los nodos internos a las rectas
            newVert = coord(intNode-nsides:intNode-1,:);
            for iMaster = 1:nsides/2
                vertA = newVert(iMaster,:);
                vertB = newVert(iMaster+1,:);
                bool = 0;
                while bool == 0
                    if norm((vertB-vertA)/div(iMaster)) > 1/divUnit
                        div(iMaster) = div(iMaster)+1;
                    else
                        bool = 1;
                    end
                end
                for intDiv = 1:div(iMaster)-1
                    sideVec = intDiv*(vertB-vertA)/div(iMaster);
                    sidePos = vertA+sideVec;
                    coord(intNode,:) = sidePos;
                    intNode = intNode+1;
                end
            end
            for iSlave = 1:nsides/2
                if iSlave == nsides/2
                    vertA = newVert(end,:);
                    vertB = newVert(1,:);
                else
                    vertA = newVert(iMaster+iSlave,:);
                    vertB = newVert(iMaster+iSlave+1,:);
                end
                for intDiv = 1:div(iSlave)-1
                    sideVec = intDiv*(vertB-vertA)/div(iSlave);
                    sidePos = vertA+sideVec;
                    coord(intNode,:) = sidePos;
                    intNode = intNode+1;
                end
            end
            div = div-1;
        end
    end
end

function connec = computeConnectivities(coord)
    connec = delaunay(coord);
end

function masterSlave = obtainMasterSlaveNodes(vert,boundary,nsides,div,dim,boundNodes)
    masterSlave = computeMasterSlaveNodes(vert,boundary,nsides,div,dim,boundNodes);
end

function plotCoordinates(coord,connec)
    s.coord = coord;
    s.connec = connec;
    m = Mesh(s);
    m.plot();
end
    
function plotVertices(vertexIndex,coord)
    plotNodes(vertexIndex,coord,'blue')
end

function plotMasterSlaveNodes(masterSlaveIndex,coord)
    masterIndex = masterSlaveIndex(:,1);
    slaveIndex  = masterSlaveIndex(:,2);
    plotNodes(masterIndex,coord,'green')
    plotNodes(slaveIndex,coord,'red')
end

function plotNodes(ind,coord,colorValue)
    b = num2str(ind);
    c = cellstr(b);
    dx = 0.01; dy = 0.01;
    x = coord(ind,1)';
    y = coord(ind,2)';
    t = text(x+dx,y+dy,c);
    set(t,'Color',colorValue)
end

function writeFEMreadingfunction(coord,connec,masterSlaveNodes,meshfilename,vertCoord)
    Data_prb = {'''TRIANGLE''','''SI''','''2D''','''Plane_Stress''','''ELASTIC''','''MICRO'''};
    
    % Initialization of the document
    fileID = fopen(meshfilename,'w');
    fprintf(fileID,'%c','%% Data file');
    fprintf(fileID,'\n\n');

    % Characteristics
    fprintf(fileID,'%c','Data_prb = {');
    fprintf(fileID,'\n');
    for k = 1:length(Data_prb)
        fprintf(fileID,'%s ;\r\n',Data_prb{k});
    end
    fprintf(fileID,'%c','};');
    fprintf(fileID,'\n\n\n');

    % Coordinates
    fprintf(fileID,'%c','coord = [');
    fprintf(fileID,'\n');
    for k = 1:size(coord,1)
        fprintf(fileID,'%g %g %g %g\r\n',k,coord(k,1),coord(k,2),0);
    end
    fprintf(fileID,'%c','];');
    fprintf(fileID,'\n\n');

    % Point loads
    fprintf(fileID,'%c','pointload = [');
    fprintf(fileID,'\n');
    fprintf(fileID,'%c','];');
    fprintf(fileID,'\n\n\n\n');

    % Connectivities
    fprintf(fileID,'%c','connec = [');
    fprintf(fileID,'\n');
    for k = 1:size(connec,1)
        fprintf(fileID,'%g %g %g %g\r\n',k,connec(k,1),connec(k,2),connec(k,3));
    end
    fprintf(fileID,'%c','];');
    fprintf(fileID,'\n\n');

    % Variable prescribed
    fprintf(fileID,'%c','%% Variable Prescribed');
    fprintf(fileID,'\n');
    fprintf(fileID,'%s \t %s \t %s','% Node','Dimension','Value');
    fprintf(fileID,'\n\n');
    fprintf(fileID,'%c','dirichlet_data = [');
    fprintf(fileID,'\n');
    for k = 1:size(vertCoord,1)
        fprintf(fileID,'%g %g %g\r\n',k,1,0);
        fprintf(fileID,'%g %g %g\r\n',k,2,0);
    end
    fprintf(fileID,'%c','];');
    fprintf(fileID,'\n\n\n');

    % Force presecribed
    fprintf(fileID,'%c','%% Force Prescribed');
    fprintf(fileID,'\n');
    fprintf(fileID,'%s \t %s \t %s','% Node','Dimension','Value');
    fprintf(fileID,'\n\n');
    fprintf(fileID,'%c','pointload_complete = [');
    fprintf(fileID,'\n');
    fprintf(fileID,'%c','];');
    fprintf(fileID,'\n\n\n');

    % Volumetric force
    fprintf(fileID,'%c','%% Volumetric Force');
    fprintf(fileID,'\n');
    fprintf(fileID,'%s \t %s \t %s','% Element','Dimension','Force_Dim');
    fprintf(fileID,'\n\n');
    fprintf(fileID,'%c','Vol_force = [');
    fprintf(fileID,'\n');
    fprintf(fileID,'%c','];');
    fprintf(fileID,'\n\n\n');

    % Group elements
    fprintf(fileID,'%c','%% Group Elements');
    fprintf(fileID,'\n');
    fprintf(fileID,'%s \t %s','% Element','Group_num');
    fprintf(fileID,'\n\n');
    fprintf(fileID,'%c','Group = [');
    fprintf(fileID,'\n');
    fprintf(fileID,'%c','];');
    fprintf(fileID,'\n\n\n');

    % Initial holes
    fprintf(fileID,'%c','%% Group Elements');
    fprintf(fileID,'\n');
    fprintf(fileID,'%c','% Elements that are considered holes initially');
    fprintf(fileID,'\n');
    fprintf(fileID,'%s','% Element');
    fprintf(fileID,'\n\n');
    fprintf(fileID,'%c','Initial_holes = [');
    fprintf(fileID,'\n');
    fprintf(fileID,'%c','];');
    fprintf(fileID,'\n\n\n');

    % Boundary elements
    fprintf(fileID,'%c','%% Boundary Elements');
    fprintf(fileID,'\n');
    fprintf(fileID,'%c','% Elements that can not be removed');
    fprintf(fileID,'\n');
    fprintf(fileID,'%s','% Element');
    fprintf(fileID,'\n\n');
    fprintf(fileID,'%c','Boundary_elements = [');
    fprintf(fileID,'\n');
    fprintf(fileID,'%c','];');
    fprintf(fileID,'\n\n\n');

    % Boundary elements
    fprintf(fileID,'%c','%% Micro Gauss Post');
    fprintf(fileID,'\n');
    fprintf(fileID,'%s','% Element');
    fprintf(fileID,'\n\n');
    fprintf(fileID,'%c','Micro_gauss_post = [');
    fprintf(fileID,'\n');
    fprintf(fileID,'%c','];');
    fprintf(fileID,'\n\n\n');

    % Master-slave nodes
    fprintf(fileID,'%c','%% Master-Slave');
    fprintf(fileID,'\n\n');
    fprintf(fileID,'%c','Master_slave = [');
    fprintf(fileID,'\n');
    for k = 1:size(masterSlaveNodes,1)
        fprintf(fileID,'%g %g\r\n',masterSlaveNodes(k,1),masterSlaveNodes(k,2));
    end
    fprintf(fileID,'%c','];');
    fprintf(fileID,'\n\n');

    % Nodes solid
    fprintf(fileID,'%c','%% Nodes Solid');
    fprintf(fileID,'\n');
    fprintf(fileID,'%c','% Nodes that must remain');
    fprintf(fileID,'\n');
    fprintf(fileID,'%s','% Nodes');
    fprintf(fileID,'\n\n');
    fprintf(fileID,'%c','% nodesolid = 1;');
    fprintf(fileID,'\n\n\n');

    % External border elements
    fprintf(fileID,'%c','%% External Border Elements');
    fprintf(fileID,'\n');
    fprintf(fileID,'%c','% Detect the elements that define the edge of the domain');
    fprintf(fileID,'\n');
    fprintf(fileID,'%s \t %s \t %s','% Element','Node(1)','Node(2)');
    fprintf(fileID,'\n\n');
    fprintf(fileID,'%c','External_border_elements = [');
    fprintf(fileID,'\n');
    fprintf(fileID,'%c','];');
    fprintf(fileID,'\n\n\n');

    % External border nodes
    fprintf(fileID,'%c','%% External Border Nodes');
    fprintf(fileID,'\n');
    fprintf(fileID,'%c','% Detect the nodes that define the edge of the domain');
    fprintf(fileID,'\n');
    fprintf(fileID,'%s \t %s \t %s','% Node');
    fprintf(fileID,'\n\n');
    fprintf(fileID,'%c','External_border_nodes = [');
    fprintf(fileID,'\n');
    fprintf(fileID,'%c','];');
    fprintf(fileID,'\n\n\n');

    % Materials
    fprintf(fileID,'%c','%% Materials');
    fprintf(fileID,'\n');
    fprintf(fileID,'%c','% Materials that have been used');
    fprintf(fileID,'\n');
    fprintf(fileID,'%s \t %s \t %s \t %s','% Material_Num','Mat_density','Young_Modulus','Poisson');
    fprintf(fileID,'\n\n');
    fprintf(fileID,'%c','Materials = [');
    fprintf(fileID,'\n');
    fprintf(fileID,'%c','];');
end



