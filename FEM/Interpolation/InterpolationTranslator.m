classdef InterpolationTranslator < handle
    
    properties (GetAccess = public, SetAccess = private)
        dim
        coord
        connec
        meshConnec
        fieldConnec
        inputBC
    end

    properties (Access = private)
        mesh
        interpolation
    end
    
    methods (Access = public)

        function obj = InterpolationTranslator(cParams)
            obj.init(cParams);
            obj.updateInputMismatch();
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
%             obj.dim           = cParams.dim;
            obj.mesh          = cParams.mesh;
            obj.interpolation = cParams.interpolation;
            obj.inputBC       = cParams.inputBC;
            obj.fieldConnec   = cParams.mesh.connec;
            obj.coord   = cParams.mesh.coord;
            obj.connec  = cParams.mesh.connec;
        end

        function updateInputMismatch(obj)
            if isequal(obj.interpolation.order, 'QUADRATIC')
                obj.updateConnectivities();
                obj.updateBoundaryConditions();
%                 obj.updateDimensions();
            end
        end

        function updateConnectivities(obj)
            s.mesh          = obj.mesh;
            s.interpolation = obj.interpolation;
            c = ConnecCoordFromInterpAndMesh(s);
            c.compute();
            obj.meshConnec  = obj.mesh.connec;
            obj.fieldConnec = c.connec;
            obj.coord       = c.coord;
            obj.connec      = c.connec;
        end

        function updateBoundaryConditions(obj)
            trimmedGlobal  = obj.fieldConnec(:,1:3);
            uniqueOldConnec = unique(obj.meshConnec, 'stable');
            uniqueNewConnec = unique(trimmedGlobal, 'stable');
            dirichletNodes = obj.inputBC.dirichlet(:,1)';
            neumannNodes  = obj.inputBC.pointload(:,1)';

            [idxD,~] = find(uniqueOldConnec == dirichletNodes);
            newDirichlet = uniqueNewConnec(idxD);
            obj.inputBC.dirichlet(:,1) = newDirichlet;
            [idxN,~] = find(uniqueOldConnec == neumannNodes);
            newNeumann = uniqueNewConnec(idxN);
            obj.inputBC.pointload(:,1) = newNeumann;
        end

%         function updateDimensions(obj)
%             obj.dim.nnodes = max(max(obj.globalConnec));
%             obj.dim.ndof   = obj.dim.nnodes * obj.dim.ndimField;
%         end

    end

end

