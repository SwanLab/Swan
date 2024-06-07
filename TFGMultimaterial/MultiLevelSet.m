classdef MultiLevelSet < handle
    
    properties (Access = public)
        designVariable
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
           % obj.createUnfittedMesh();
            obj.create();
        end

        function update(obj)
            levelSet = obj.designVariable.value; 
            levelSet = levelSet./normL2(obj.unitM,levelSet ); % level-set function nomalization  
            obj.designVariable.update(levelSet)
        end

        function updateOld(obj)
            obj.funOld = obj.lsFun; % .copy();
        end

        function charFun = obtainDomainFunction(obj)
            uMesh = obj.getUnfittedMesh();
            m = obj.mesh;
            charFun = CharacteristicFunction.create(m);
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
            s.order                = 'P1';
            s.mesh                 = obj.mesh;
            s.fun                  = LagrangianFunction(s);
            ls1                    = DesignVariable.create(s);
            obj.designVariable{1} = ls1;

            s.fValues              = obj.lsFun{2};
            s.fun                  = LagrangianFunction(s);
            ls1                    = DesignVariable.create(s);
            obj.designVariable{2} = ls1;

            s.fValues              = obj.lsFun{3};
            s.fun                  = LagrangianFunction(s);
            ls1                    = DesignVariable.create(s);
            obj.designVariable{3} = ls1;
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