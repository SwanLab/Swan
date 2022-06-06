classdef InterpolationTranslator < handle
    
    properties (GetAccess = public, SetAccess = private)
        coord
        connec
        meshConnec
        fieldConnec
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

        function newBC = updateBoundaryConditions(obj, inputBC)
            newBC = inputBC;
            if isequal(obj.interpolation.order, 'QUADRATIC')
                trimmedGlobal  = obj.fieldConnec(:,1:3);
                uniqueOldConnec = unique(obj.meshConnec, 'stable');
                uniqueNewConnec = unique(trimmedGlobal, 'stable');
                dirichletNodes = inputBC.dirichlet(:,1)';
                [idxD,~] = find(uniqueOldConnec == dirichletNodes);
                newDirichlet = uniqueNewConnec(idxD);
                inputBC.dirichlet(:,1) = newDirichlet;
                if ~isempty(inputBC.pointload)
                    neumannNodes  = inputBC.pointload(:,1)';
                    [idxN,~] = find(uniqueOldConnec == neumannNodes);
                    newNeumann = uniqueNewConnec(idxN);
                    inputBC.pointload(:,1) = newNeumann;
                end
                newBC = inputBC;
            end
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh          = cParams.mesh;
            obj.interpolation = cParams.interpolation;
            obj.fieldConnec   = cParams.mesh.connec;
            obj.coord   = cParams.mesh.coord;
            obj.connec  = cParams.mesh.connec;
        end

        function updateInputMismatch(obj)
            if isequal(obj.interpolation.order, 'QUADRATIC')
                obj.updateConnectivities();
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

    end

end

