classdef MultiLevelSet < handle
    
    properties (Access = public)
        levelSets
        fun
    end

    properties (Access = private)
        type
        lsFun
        funOld
        plotting
        mesh
        unitM
        unfittedMesh
    end

    methods (Access = public)

        function obj = MultiLevelSet(cParams)
            obj.init(cParams);
            obj.create();
            obj.plotting = cParams.plotting;
        end

        function update(obj,value)
            for i = 1:size(value,2)
                obj.levelSets{1,i}.update(value(:,i));
            end
            obj.fun.fValues = [];
            obj.fun.fValues = [obj.fun.fValues;obj.levelSets{1}.fun.fValues];
            obj.fun.fValues = [obj.fun.fValues;obj.levelSets{2}.fun.fValues];
            obj.fun.fValues = [obj.fun.fValues;obj.levelSets{3}.fun.fValues];
        end

        function plot(obj)
            if obj.plotting
                charfun = obj.obtainDomainFunction();
                [fi,~]  = charfun.computeAtNodesAndElements();
                p          = obj.mesh.coord';
                t          = obj.mesh.connec';
                t(4,:)     = 1;

                figure(2)
                multimat_plot(p,t,fi.fValues);
                drawnow
            end
        end

        function updateOld(obj)
            obj.funOld = obj.lsFun; % .copy();
            for i = 1:length(obj.levelSets)
                obj.levelSets{i}.updateOld();
            end
        end

        function norm = computeL2normIncrement(obj)
            norm = 0;
            for i = 1:length(obj.levelSets)
                ni = obj.levelSets{i}.computeL2normIncrement();
                norm = norm + ni^2;
            end
            norm = sqrt(norm);
        end

        function charFun = obtainDomainFunction(obj)
            s.designVariable = obj;
            s.mesh           = obj.mesh; 
            charFun          = MultiMaterialCharacteristicFunction(s);
        end

        function m = getUnfittedMesh(obj)
            m = obj.unfittedMesh;
        end
    end

   

    methods (Access = private)
        function init(obj,cParams)
            obj.type = cParams.type;
            obj.mesh = cParams.mesh;
            obj.lsFun  = cParams.lsFun;
            obj.unitM = cParams.unitM;
        end

        function create(obj)
            s.type                 = 'LevelSet';
            s.plotting             = false;
            s.fValues              = obj.lsFun{1};
            s.mesh                 = obj.mesh;
            s.ndimf                = 1;
            s.fHandle              = obj.lsFun{1};
            s.fun                  = AnalyticalFunction(s);
            s.fun                  = s.fun.project('P1');
            ls1                    = DesignVariable.create(s);
            obj.levelSets{1} = ls1;

            s.fHandle              = obj.lsFun{2};
            s.fun                  = AnalyticalFunction(s);
            s.fun                  = s.fun.project('P1');
            ls2                    = DesignVariable.create(s);
            obj.levelSets{2} = ls2;

            s.fHandle              = obj.lsFun{3};
            s.fun                  = AnalyticalFunction(s);
            s.fun                  = s.fun.project('P1');
            ls3                    = DesignVariable.create(s);
            obj.levelSets{3} = ls3;

            obj.fun.fValues = [];
            obj.fun.fValues = [obj.fun.fValues;obj.levelSets{1}.fun.fValues];
            obj.fun.fValues = [obj.fun.fValues;obj.levelSets{2}.fun.fValues];
            obj.fun.fValues = [obj.fun.fValues;obj.levelSets{3}.fun.fValues];
        end

        function createUnfittedMesh(obj)
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            cParams = SettingsMeshUnfitted(s);
            obj.unfittedMesh = UnfittedMesh(cParams);
            obj.updateUnfittedMesh();
        end
        
        function updateUnfittedMesh(obj)
            obj.unfittedMesh.compute(obj.lsFun); %.fValues);
        end
    end

end