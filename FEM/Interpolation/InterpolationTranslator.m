classdef InterpolationTranslator < handle
    
    properties (GetAccess = public, SetAccess = private)
        dim
        coord
        linealConnec
        globalConnec
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
            obj.dim           = cParams.dim;
            obj.mesh          = cParams.mesh;
            obj.interpolation = cParams.interpolation;
            obj.inputBC       = cParams.inputBC;
            obj.globalConnec  = cParams.mesh.connec;
        end

        function updateInputMismatch(obj)
            if isequal(obj.interpolation.order, 'QUADRATIC')
                obj.updateConnectivities();
                obj.updateBoundaryConditions();
                obj.updateDimensions();
            end
        end

        function updateConnectivities(obj)
            s.mesh          = obj.mesh;
            s.interpolation = obj.interpolation;
            c = ConnecCoordFromInterpAndMesh(s);
            c.compute();
            obj.linealConnec = obj.mesh.connec;
            obj.globalConnec = c.connec;
            obj.coord        = c.coord;
        end

        function updateBoundaryConditions(obj)
            trimmedGlobal  = obj.globalConnec(:,1:3);
            uniqueOldConnec = unique(obj.linealConnec, 'stable');
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

        function updateDimensions(obj)
            obj.dim.nnodes = max(max(obj.globalConnec));
            obj.dim.ndof   = obj.dim.nnodes * obj.dim.ndimField;
        end

    end

end

