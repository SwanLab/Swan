classdef NodesConnector < handle
    
    properties (Access = private)
        nodes
        coord
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
            if (obj.nodes.vert == 6) %|| (obj.theta(2)~=90)
                obj.deleteExtraElements();
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nodes = cParams.nodes;
            obj.coord = cParams.coord;
            obj.theta = cParams.theta;
        end
        
        function deleteExtraElements(obj)
            cont = 1;
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
        
    end
    
end