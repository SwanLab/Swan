classdef DimensionScalar < handle
    
    properties (GetAccess = public, SetAccess = private)
        name
        nelem
        npnod
        nnode
        ndimField
        ndofPerElement
        ndof

        nstre %nvoigt
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
%             obj.mesh           = cParams.mesh;
            obj.ndimField      = 1; % by definition
            obj.name           = cParams.name;
            obj.nnode          = cParams.mesh.interpolation.nnode;
            obj.ndofPerElement = cParams.mesh.interpolation.nnode;
            obj.ndof           = cParams.mesh.npnod;

            obj.nelem          = cParams.mesh.nelem;
            obj.npnod          = cParams.mesh.npnod;

            obj.nstre = 2; % by definition. does it make sense?
        end

    end

end

