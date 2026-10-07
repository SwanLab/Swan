classdef NodesConnector < handle
    
    properties (Access = private)
        nodes
        coord
        div
        latticeVectors
    end
    
    properties (Access = public)
        connec
    end
    
    methods (Access = public)

        function obj = NodesConnector(cParams)
            obj.init(cParams);
        end

        function computeConnections(obj)

            if obj.nodes.vert == 6

                % Keep the original implementation for hexagonal cells
                obj.connec = delaunay(obj.coord);
                obj.deleteExtraElementsCaseA();

            elseif obj.nodes.vert == 4

                % Structured connectivity for parallelogram cells
                obj.computeStructuredParallelogramConnectivity();

            else
                error('Unsupported number of vertices.')
            end

        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nodes = cParams.nodes;
            obj.coord = cParams.coord;
            obj.div   = cParams.div;
            if isfield(cParams, 'latticeVectors')
                obj.latticeVectors = cParams.latticeVectors;
            end
        end

        function deleteExtraElementsCaseA(obj)
            cont = 1;
            rowsToDelete = [];
            for i = 1:size(obj.connec,1)
                found = zeros(1,3);
                elementNodes = obj.connec(i,:);
                for j = 1:size(obj.connec,2)
                    if elementNodes(j) < obj.nodes.bound+1
                        found(j) = found(j)+1;
                    end
                end
                if sum(found) == 3
                    rowsToDelete(cont,1) = i;
                    cont = cont+1;
                end
            end
            obj.connec(rowsToDelete,:) = [];
        end
        function computeStructuredParallelogramConnectivity(obj)
            n1 = obj.div(1);
            n2 = obj.div(2);
            nVert = obj.nodes.vert;
            nBound = obj.nodes.bound;
            startL1 = nVert + 1;
            startL2 = startL1 + (n1 - 1);
            startL3 = startL2 + (n2 - 1);
            startL4 = startL3 + (n1 - 1);
            nodeId = zeros(n2+1,n1+1);
            nodeId(1,1)       = 1;   
            nodeId(1,n1+1)    = 2;  
            nodeId(n2+1,n1+1) = 3;   
            nodeId(n2+1,1)    = 4; 
            for i = 1:n1-1
                nodeId(1,i+1) = startL1 + (i-1);
            end
            for j = 1:n2-1
                nodeId(j+1,n1+1) = startL2 + (j-1);
            end
            for i = 1:n1-1

                nodeId(n2+1,i+1) = startL3 + (n1-i-1);
            end
            for j = 1:n2-1
                nodeId(j+1,1) = startL4 + (n2-j-1);
            end
            for j = 1:n2-1
                for i = 1:n1-1
                    nodeId(j+1,i+1) = nBound + (j-1)*(n1-1) + i;
                end
            end
            nElem = 2*n1*n2;
            connec = zeros(nElem,3);
            e = 1;
            for j = 1:n2
                for i = 1:n1
                    n00 = nodeId(j,  i);
                    n10 = nodeId(j,  i+1);
                    n01 = nodeId(j+1,i);
                    n11 = nodeId(j+1,i+1);
                    connec(e,:) = [n00 n10 n11];
                    e = e + 1;
                    connec(e,:) = [n00 n11 n01];
                    e = e + 1;
                end
            end
            obj.connec = connec;
        end        
        function deleteExtraElementsCaseB(obj)
            cont = 1;
            rowsToDelete = [];
            for i = 1:size(obj.connec,1)
                found = zeros(1,3);
                noDeletionA = zeros(1,3);
                noDeletionB = zeros(1,3);
                noDeletionC = zeros(1,3);
                noDeletionD = zeros(1,3);
                elementNodes = obj.connec(i,:);
                for j = 1:size(obj.connec,2)
                    if (elementNodes(j) < obj.nodes.bound+1)
                        found(j) = found(j)+1;
                        if (elementNodes(j) == 1) || (elementNodes(j) == obj.nodes.vert+1) || (elementNodes(j) == obj.nodes.bound)  
                            noDeletionA(j) = noDeletionA(j)+1;
                        end
                        if (elementNodes(j) == 2) || (elementNodes(j) == obj.nodes.vert+obj.div(1)-1) || (elementNodes(j) == obj.nodes.vert+obj.div(1))
                            noDeletionB(j) = noDeletionB(j)+1;
                        end
                        if (elementNodes(j) == 3) || (elementNodes(j) == obj.nodes.vert+obj.div(1)+obj.div(2)-2) || (elementNodes(j) == obj.nodes.vert+obj.div(1)+obj.div(2)-1)
                            noDeletionC(j) = noDeletionC(j)+1;
                        end
                        if (elementNodes(j) == 4) || (elementNodes(j) == obj.nodes.vert+2*obj.div(1)+obj.div(2)-3) || (elementNodes(j) == obj.nodes.vert+2*obj.div(1)+obj.div(2)-2)
                            noDeletionD(j) = noDeletionD(j)+1;
                        end
                    end
                end
                if (sum(found) == 3) && (sum(noDeletionA) ~= 3) && (sum(noDeletionB) ~= 3) && (sum(noDeletionC) ~= 3) && (sum(noDeletionD) ~= 3)
                    rowsToDelete(cont,1) = i;
                    cont = cont+1;
                end
            end
            obj.connec(rowsToDelete,:) = [];
        end
        
    end
    
end