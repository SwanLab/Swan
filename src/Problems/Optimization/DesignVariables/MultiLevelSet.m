classdef MultiLevelSet < handle
    
    properties (Access = public)
        levelSets
        fun
    end

    properties (Access = private)
        plotting
        mesh
    end

    methods (Access = public)
        function obj = MultiLevelSet(cParams)
            obj.init(cParams);
            obj.create(cParams);
        end

        function update(obj,value)
            obj.fun.fValues = value;
            nLS             = length(obj.levelSets);
            value           = reshape(value,[],nLS);
            for i = 1:size(value,2)
                obj.levelSets{1,i}.update(value(:,i));
            end
        end

        function plot(obj)
            if obj.plotting
            % WE HAVE TO IMPLEMENT OUR VERSION
            end
        end

        function updateOld(obj)
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

        function charFuns = obtainDomainFunction(obj)
            nLS = length(obj.levelSets);
            charFuns = cell(nLS,1);
            for i = 1:nLS
                charFuns(i) = obj.levelSets{i}.obtainDomainFunction();
            end
        end

        function charFun = obtainGlobalDomainFunction(obj)
            s.designVariable = obj;
            s.mesh           = obj.mesh; 
            multiCharFun     = MultiMaterialCharacteristicFunction(s);
            charFun          = multiCharFun.compute();
        end

        function charFun = obtainDomainFunctionDerivatives(obj)
            s.designVariable = obj;
            s.mesh           = obj.mesh;
            multiCharFun     = MultiMaterialCharacteristicFunction(s);
            charFun          = multiCharFun.computeDerivatives();
        end

        function ls = obtainVariableInCell(obj)
            ls = obj.levelSets;
        end
    end

   

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh     = cParams.mesh;
            obj.plotting = cParams.plotting;
        end

        function create(obj,s)
            for i = 1:length(s.lsFun)
                obj.levelSets{i} = obj.createLevelSet(s.lsFun{i});
            end
            obj.createFun();
        end

        function ls = createLevelSet(obj,lsFun)
            s.fun      = lsFun;
            s.mesh     = obj.mesh;
            s.type     = 'LevelSet';
            s.plotting = false;
            ls         = DesignVariable.create(s);
        end

        function createFun(obj)
            obj.fun.fValues = [];
            for i = 1:length(obj.levelSets)
                obj.fun.fValues = [obj.fun.fValues;obj.levelSets{i}.fun.fValues];
            end
        end
    end
end