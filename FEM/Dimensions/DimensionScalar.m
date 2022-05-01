classdef DimensionScalar < handle
    
    properties (GetAccess = public, SetAccess = private)
        name
        npnod
        nnodeElem
        ndimField
        ndofPerElement
        ndof
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
            obj.ndimField      = 1; % by definition
            obj.name           = cParams.name;
            obj.nnodeElem      = cParams.mesh.interpolation.nnode;
            obj.ndofPerElement = cParams.mesh.interpolation.nnode;
            obj.ndof           = cParams.mesh.npnod;

            obj.npnod          = cParams.mesh.npnod;
        end

    end

end

