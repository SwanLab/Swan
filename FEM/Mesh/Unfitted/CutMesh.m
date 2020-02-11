classdef CutMesh < Mesh
    
    properties (GetAccess = public, SetAccess = private)
         subcellIsoCoords
         cellContainingSubcell
         globalConnec
    end
    
    properties (Access = private)
        backgroundMesh
    end
        
    methods (Access = public)
        
        function obj = CutMesh(cParams)
            obj.coord  = cParams.coord;
            obj.connec = cParams.connec;
            obj.type   = cParams.type;
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.subcellIsoCoords = cParams.subcellIsoCoords;
            obj.cellContainingSubcell = cParams.cellContainingSubcell;
            obj.computeDescriptorParams();
            obj.createInterpolation();
            obj.computeElementCoordinates();    
            obj.computeGlobalConnec();
        end
        
    end      
    
    methods (Access = protected)
        
        function computeEmbeddingDim(obj)
            switch obj.type
                case 'BOUNDARY'
                    obj.embeddedDim = obj.ndim - 1;
                case {'INTERIOR','COMPOSITE'}
                    obj.embeddedDim = obj.ndim;
                otherwise
                    error('EmbeddedDim not defined')
            end
        end
        
    end 
    
    methods (Access = private)
        
        function computeGlobalConnec(obj)
            nnode = obj.backgroundMesh.nnode;
            nelem = obj.nelem;
            obj.globalConnec = zeros(nelem,nnode);
            for ielem = 1:nelem
                icell  = obj.cellContainingSubcell(ielem);
                nodes  = obj.backgroundMesh.connec(icell,:);
                obj.globalConnec(ielem,:) = nodes;
            end            
            
        end
    end
    
end