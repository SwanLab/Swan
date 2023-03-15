classdef LevelSet < DesignVariable
    
    properties (Access = private)
        creatorSettings
        unfittedMesh
    end
    
    methods (Access = public)
        
        function obj = LevelSet(cParams)
            obj.nVariables = 1;
            obj.init(cParams);
            obj.creatorSettings      = cParams.creatorSettings;
            obj.creatorSettings.type = cParams.initialCase;
            obj.createValue();
            obj.createUnfittedMesh();
        end
        
        function update(obj,value)
            obj.value = value;
            if ~isempty(obj.isFixed)
                obj.value(obj.isFixed.nodes) = obj.isFixed.values;
            end
            obj.updateUnfittedMesh();
        end
        
        function m = getUnfittedMesh(obj)
            m = obj.unfittedMesh;
        end
        
        function v = getVariablesToPlot(obj)
            v{1} = obj.value;
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
        
        function [fun, funNames] = getFunsToPlot(obj)
            aa.mesh = obj.mesh.meshes{1};
            aa.fValues = obj.value;
            valFun = P1Function(aa);
            
%             bb.mesh = obj.mesh.meshes{1};
%             bb.fValues = obj.alpha';
%             alphaFun = P0Function(bb);
% 
%             fun = {valFun, alphaFun};
%             funNames = {'value', 'alpha'};
            fun = {valFun};
            funNames = {'value'};
            cc.mesh     = obj.mesh.meshes{1};
            cc.filename = 'levelsetdesignvar';
            cc.fun      = fun;
            cc.funNames = funNames;
%             pvPst = ParaviewPostprocessor(cc);
%             pvPst.print();
%             fp = FunctionPrinter(cc);
%             fp.print();
        end
    end
    
    methods (Access = private)
        
        function createValue(obj)
            s         = obj.creatorSettings;
            s.ndim    = obj.mesh.ndim;
            s.coord   = obj.mesh.coord;
            lsCreator = LevelSetCreator.create(s);
            lsValue   = lsCreator.getValue();
            obj.value = lsValue;
        end
        
        function createUnfittedMesh(obj)
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            cParams = SettingsMeshUnfitted(s);
            obj.unfittedMesh = UnfittedMesh(cParams);
            obj.updateUnfittedMesh();
        end
        
        function updateUnfittedMesh(obj)
            obj.unfittedMesh.compute(obj.value);
        end
        
    end
    
end