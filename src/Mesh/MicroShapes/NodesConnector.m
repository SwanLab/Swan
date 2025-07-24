classdef NodesConnector < handle
    
    properties (Access = private)
        nodes
        coord
        div
        theta
    end
    
    properties (Access = public)
        connec
    end
    
    methods (Access = public)
        
        function obj = NodesConnector(cParams)
            obj.init(cParams);
        end
        
        function computeConnections(obj)
            obj.connec = delaunay(obj.coord);
            if obj.nodes.vert == 6
                obj.deleteExtraElementsCaseA();
            elseif obj.theta(2) ~= 90
                obj.deleteExtraElementsCaseB();
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nodes = cParams.nodes;
            obj.coord = cParams.coord;
            obj.theta = cParams.theta;
            obj.div = cParams.div;
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