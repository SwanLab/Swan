classdef SubcellsMesher_Boundary_3D < SubcellsMesher_Boundary
    
    properties (Access = private)
        subcells_connec
        boundary_subcells_connec
        nBoundSubcells
        nIntNodes
        nExtNodes
    end
    
    methods (Access = public)
        
        function obj = SubcellsMesher_Boundary_3D(cParams)
           obj.init(cParams); 
        end
        
    end        
        
    methods (Access = protected)
        
        function computeFacetsConnectivities(obj)
            obj.computeAllPossibleSubcellsInCell();
            obj.getBoundarySubcells();
            obj.removeCellNodesFromBoundarySubcells();
        end
        
    end
    
    methods (Access = private)
        
        function computeAllPossibleSubcellsInCell(obj)
            obj.subcells_connec = obj.computeDelaunay(obj.interior_coord_iso);
        end
        
        function getBoundarySubcells(obj)
            obj.classifyNodes();
            intCond = obj.nIntNodes == 1;
            extCond = obj.nExtNodes == 0;
            obj.boundary_subcells_connec = obj.subcells_connec(intCond & extCond,:);
            obj.nBoundSubcells = size(obj.boundary_subcells_connec,1);
        end
        
        function removeCellNodesFromBoundarySubcells(obj)
            obj.allocateMemoryConnec();
            for i = 1:obj.nBoundSubcells
                cutPoints = obj.boundary_subcells_connec(i,:)>obj.nCellNodes;
                obj.connec(i,:) = obj.boundary_subcells_connec(i,cutPoints);
            end
            obj.connec = obj.connec - obj.nCellNodes;
        end
        
        function classifyNodes(obj)
            phi = obj.cell_levelSet;
            
            intNodes = find(phi<=0);
            extNodes = find(phi>0);
            
            obj.nIntNodes = obj.countInputNodesPerCell(obj.subcells_connec,intNodes); %#ok<FNDSB>
            obj.nExtNodes = obj.countInputNodesPerCell(obj.subcells_connec,extNodes); %#ok<FNDSB>
        end
        
        function allocateMemoryConnec(obj)
            %nnode = size(obj.connec,2);
            nnode = 3;
            obj.connec = zeros([obj.nBoundSubcells,nnode]);
        end
        
    end
    
    methods (Access = private, Static)
        
        function n = countInputNodesPerCell(connec,nodes)
            n = zeros(size(connec,1),1);
            for inode = 1:length(nodes)
                match = false(size(connec,1),1);
                for iconnec = 1:size(connec,2)
                    match = match | connec(:,iconnec) == nodes(inode);
                end
                n = n + match ;
            end
        end
        
    end
    
end

