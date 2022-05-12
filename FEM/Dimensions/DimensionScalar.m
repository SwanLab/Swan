classdef DimensionScalar < handle
    
    properties (GetAccess = public, SetAccess = private)
        name
        nnodes
        nnodeElem
        ndimf
        ndofsElem
        ndofs
    end
    
    properties (Access = private)
%         mesh
    end
    
    methods (Access = public)

        function obj = DimensionScalar(cParams)
            obj.init(cParams);
        end
        
    end

    methods (Access = private)

        function init(obj, cParams)
            obj.ndimf          = 1; % by definition
            obj.name           = cParams.name;
            obj.nnodeElem      = cParams.mesh.interpolation.nnode;
            obj.ndofsElem      = cParams.mesh.interpolation.nnode;
            obj.ndofs          = cParams.mesh.nnodes;
            obj.nnodes         = cParams.mesh.nnodes;
        end

    end

end

