function masterSlave = computeMasterSlaveNodes(vert,bound,nsides,div,dim,boundNodes)
    % Init of the initial data (function)
    [normalVec,ortoNodes] = computeBoundaryMeshes(vert,bound,nsides,div,dim);
    pairs = computePairOfMeshes(normalVec,dim,nsides);
    masterSlave = computeMasterAndSlaves(ortoNodes,pairs,dim,div,nsides,boundNodes);
end


function [normalVec,ortoNodes] = computeBoundaryMeshes(vert,bound,nsides,div,dim)
    normalVec = zeros(nsides,2);
    for iVert = 1:nsides
        % función de obtención de vertices
        vertexA = vert(iVert,:);
        if iVert == nsides
            vertexB = vert(1,:);
        else
            vertexB = vert(iVert+1,:);
        end
        % función de obtención del vector unitario
        vec = vertexB - vertexA;
        mod = norm(vec);
        vecU = vec/mod;
        % función de obtención de la normal (siempre saliente al lado)
        normalVec(iVert,1) = normalVec(iVert,1)+vecU(1,2);
        normalVec(iVert,2) = normalVec(iVert,2)-vecU(1,1);
    end

    % Función que calcula todos los vectores posibles hasta el verticeA
    %nIntNodes = nsides*(div-1);
    nIntNodes = sum(2*(div-1));
    intNodes = bound(nsides+1:end,:);
    coord2A = zeros(nIntNodes,dim,nsides);
    for iVert = 1:nsides
        vertexA = vert(iVert,:);
        for iVec = 1:nIntNodes
            vertexC = intNodes(iVec,:);
            % función de obtención del vector unitario
            vec = vertexC - vertexA;
            mod = norm(vec);
            vecU = vec/mod;
            coord2A(iVec,:,iVert) = coord2A(iVec,:,iVert)+vecU;
        end
    end

    % find nodes with coord2A ortogonal to normal --> ortoNodes
    tol = 10e-6;
    ortoNodes = zeros(max(div)-1,nsides);
    for iVert = 1:nsides
        iPos = 1;
        vectorA = normalVec(iVert,:);
        for iVec = 1:nIntNodes
            vectorB = coord2A(iVec,:,iVert);
            if abs(dot(vectorA,vectorB)) < tol
                nodeNumber = nsides+iVec;
                ortoNodes(iPos,iVert) = ortoNodes(iPos,iVert)+nodeNumber;
                iPos = iPos+1;
            end
        end
        % Funcionalidad que voltea los nodos de los lados opuestos a los master
        if iVert > nsides/2
            ortoNodes(:,iVert) = sort(ortoNodes(:,iVert),'descend');
        end
    end
end

function pairs = computePairOfMeshes(normalVec,dim,nsides)
    tol = 10e-6;
    pairs = zeros(nsides/2,dim);
    irow = 1;
    for iSrch = 1:nsides/2
        jSrch = iSrch+1;
        vectorA = normalVec(iSrch,:);
        found = 0;
        while found == 0
            vectorB = normalVec(jSrch,:);
            cosAngle = abs(dot(vectorA,vectorB));
            if abs(cosAngle-1) < tol
                pairs(irow,1) = pairs(irow,1)+iSrch;
                pairs(irow,2) = pairs(irow,2)+jSrch;
                irow = irow+1;
                found = 1;
            else
                jSrch = jSrch+1;
            end
        end
    end
end

function masterSlave = computeMasterAndSlaves(ortoNodes,pairs,dim,div,nsides,boundNodes)
    masterSlave = zeros((boundNodes-nsides)/2,dim);
    cont = 1;
    for iPair = 1:nsides/2
        lineA = pairs(iPair,1);
        lineB = pairs(iPair,2);
        for iDiv = 1:div(iPair)-1
            masterSlave(cont,1) = masterSlave(cont,1)+ortoNodes(iDiv,lineA);
            masterSlave(cont,2) = masterSlave(cont,2)+ortoNodes(iDiv,lineB);
            cont = cont+1;
        end
    end
end
