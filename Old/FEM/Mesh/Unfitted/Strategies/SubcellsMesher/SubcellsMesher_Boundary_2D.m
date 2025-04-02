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
            type = 'QUAD';
            
            int = Interpolation.create(type,'LINEAR');
            xG = coordIso(5:end,:)';[0 0]';
            shapes = int.computeShapeFunctions(xG);
            nnode = size(shapes,1);%obj.mesh.meshBackground.nnode;
            lsNodes = levelSet(1:nnode);
            levelSetNew = shapes'*lsNodes;
            levelSetNew = zeros(4,1);
            
            levelSet = [levelSet;levelSetNew];
%             
%             m.meshBackground = obj.mesh.meshBackground;
%             m.levelSet            = m.levelSet_background;            
%             m.backgroundCutCells  = m.backgroundCutCells;
%             m.backgroundGeomInterpolation = m.backgroundGeomInterpolation;             
            
            %coordIso = [obj.interior_coord_iso(1:4,:);xG']; 
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
                connecV(idel,posCell) = connecDel(pos,posCell) - nnode;
            end
            
            lsElem(:,1) = levelSet(connecDel(:,1));
            lsElem(:,2) = levelSet(connecDel(:,2));
            lsElem(:,3) = levelSet(connecDel(:,3)); 
            
            isEmpty = any(lsElem>0,2);

            obj.connec = connecDel(~isEmpty,:);
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

