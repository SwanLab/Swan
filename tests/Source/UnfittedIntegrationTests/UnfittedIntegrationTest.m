classdef UnfittedIntegrationTest < handle
    
    properties (Access = private)
        mesh
        meshType
        levelSet
        analyticalValue
    end

    properties (Access = private)
        unfittedMesh
        varAdim
    end
    
    methods (Access = public)
        
        function obj = UnfittedIntegrationTest(cParams)
            obj.init(cParams);
            obj.integrateSurface();
        end
        
        function error = computeError(obj)
            error = abs(obj.varAdim - 1);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh            = cParams.mesh;
            obj.meshType        = cParams.meshType;
            obj.levelSet        = cParams.levelSet;
            obj.analyticalValue = cParams.analyticalValue;
        end
        
        function integrateSurface(obj)
            obj.createMesh();
            geomVar     = obj.computeGeometricalVariable();
            obj.varAdim = geomVar/obj.analyticalValue;
        end

        function createMesh(obj)
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            obj.unfittedMesh = UnfittedMesh(s);
            obj.unfittedMesh.compute(obj.levelSet);
        end

        function totalIntegral = computeGeometricalVariable(obj)
            switch obj.meshType
                case 'INTERIOR'
                    totalIntegral = obj.unfittedMesh.computeMass();
                case 'BOUNDARY'
                    totalIntegral = obj.unfittedMesh.computePerimeter();
            end
        end

    end

end