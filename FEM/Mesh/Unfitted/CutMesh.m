classdef CutMesh < Mesh
    
    properties (GetAccess = public, SetAccess = private)
         subcellIsoCoords
         cellContainingSubcell
    end
    
    properties (Access = private)
        backgroundMesh
        type
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
        end
        
        function b = getBackgroundMesh(obj)
            b = obj.backgroundMesh;
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
    
end