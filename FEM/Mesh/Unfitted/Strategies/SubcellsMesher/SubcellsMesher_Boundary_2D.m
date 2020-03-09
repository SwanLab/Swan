classdef SubcellsMesher_Boundary_2D < SubcellsMesher_Boundary
    
    methods (Access = public)
        
        function obj = SubcellsMesher_Boundary_2D(cParams)
           obj.init(cParams); 
        end
        
    end        
    
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
            levelSet = obj.cell_levelSet;
            coordIso = obj.interior_coord_iso;
            coordIso(end,:) = [0 0];
            int = Interpolation.create(obj.mesh.meshBackground,'LINEAR');
            xG = [0 0]';
            int.computeShapeDeriv(xG)
            shapes = int.shape;
            nnode = obj.mesh.meshBackground.nnode;
            lsNodes = levelSet(1:nnode);
            levelSet(end+1) = shapes'*lsNodes;
            
%             
%             m.meshBackground = obj.mesh.meshBackground;
%             m.levelSet            = m.levelSet_background;            
%             m.backgroundCutCells  = m.backgroundCutCells;
%             m.backgroundGeomInterpolation = m.backgroundGeomInterpolation;             
            
            
           % levelSet(end,1) = 
            connecDel = obj.computeDelaunay(coordIso);
            nnode = size(connecDel,2);
            
            positiveCellNodes  = find(levelSet > 0);
            nPositiveCellNodes = length(positiveCellNodes);

            connecV = zeros(nPositiveCellNodes,nnode);
            for idel = 1:nPositiveCellNodes
                [connec_positive, ~] = find(connecDel == positiveCellNodes(idel));
                nodePos = connecDel(connec_positive(end),:);
                posCell = nodePos ~= positiveCellNodes(idel);
                pos = connec_positive(end);
                connecV(idel,:) = connecDel(pos,posCell) - nnode;
            end
            obj.connec = connecV;
        end
        
        
%         
%         
%         
%         function computeCutPoints()
%             m.meshBackground = ;
%             m.levelSet            = m.levelSet_background;            
%             m.backgroundCutCells  = m.backgroundCutCells;
%             m.backgroundGeomInterpolation = m.backgroundGeomInterpolation;            
%             
%             
%         end
%         
            
            
        
        
        
        
        
        
    end
    
end

