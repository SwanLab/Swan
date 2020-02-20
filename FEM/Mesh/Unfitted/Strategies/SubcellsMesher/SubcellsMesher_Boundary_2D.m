classdef SubcellsMesher_Boundary_2D < SubcellsMesher_Boundary
    
    methods (Access = protected)
        
        function computeFacetsConnectivities(obj)
            npnod =  size(obj.coord_iso,1);
            switch npnod
                case {2,1}
                    obj.connectPairOfPoints();
                case 4
                    obj.computeConnectivitiesSolvingAmbiguity();
                otherwise
                    error('Case not considered.')
            end
            
        end
        
    end
    
    methods (Access = private)
        
        function connectPairOfPoints(obj)
            obj.connec = [1 2];
        end
        
        function computeConnectivitiesSolvingAmbiguity(obj)
            del_connec = obj.computeDelaunay(obj.interior_coord_iso);
            nnode = size(del_connec,2);
            
            positiveCellNodes = find(obj.cell_levelSet > 0);
            nPositiveCellNodes = length(positiveCellNodes);

            obj.connec = zeros(nPositiveCellNodes,nnode);
            for idel = 1:nPositiveCellNodes
                [connec_positive, ~] = find(del_connec == positiveCellNodes(idel));
                obj.connec(idel,:) = del_connec(connec_positive(end),del_connec(connec_positive(end),:)~=positiveCellNodes(idel)) - nnode;
            end
        end
    end
    
end

