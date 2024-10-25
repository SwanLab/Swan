classdef LevelSet < DesignVariable
    
    properties (Access = private)
        unfittedMesh
        plotting
        plotter
    end
    
    methods (Access = public)
        
        function obj = LevelSet(cParams)
            obj.nVariables = 1;
            obj.init(cParams);
            obj.createUnfittedMesh();
            obj.createPlotter(cParams);
        end

        function update(obj,value)
            if ~isempty(obj.isFixed)
                value(obj.isFixed.nodes) = obj.isFixed.values;
            end
            s.mesh    = obj.mesh;
            s.fValues = value;
            s.order   = 'P1';
            obj.fun   = LagrangianFunction(s);
            obj.updateUnfittedMesh();
        end

        function charFun = obtainDomainFunction(obj)
            uMesh = obj.getUnfittedMesh();
            charFun = CharacteristicFunction.create(uMesh);
        end

        function m = getUnfittedMesh(obj)
            m = obj.unfittedMesh;
        end
        
        function Vf = computeVolumeFraction(obj)
            q = Quadrature.set(obj.unfittedMesh.backgroundMesh.type);
            q.computeQuadrature('CONSTANT');
            bM = obj.unfittedMesh.backgroundMesh;
            dv = obj.unfittedMesh.computeDvolume(q);
            dVT = bM.computeDvolume(q)';
            vf = dv./dVT;
            Vf(1,1,:) = vf;
        end

        function plot(obj)
            if obj.plotting
                obj.plotter.plot();
            end
        end

        function ls = copy(obj)
            s.fun      = obj.fun;
            s.mesh     = obj.mesh;
            s.type     = 'LevelSet';
            s.plotting = false;
            ls         = DesignVariable.create(s);
        end

        function ls = obtainVariableInCell(obj)
            ls{1} = obj;
        end
    end

    methods (Access = private)

        function createUnfittedMesh(obj)
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            cParams = SettingsMeshUnfitted(s);
            obj.unfittedMesh = UnfittedMesh(cParams);
            obj.updateUnfittedMesh();
        end
        
        function updateUnfittedMesh(obj)
            obj.unfittedMesh.compute(obj.fun.fValues);
        end

        function createPlotter(obj,cParams)
            obj.plotting = cParams.plotting;
            if obj.plotting
                obj.plotter  = Plotter.create(obj);
            end
        end

    end

end