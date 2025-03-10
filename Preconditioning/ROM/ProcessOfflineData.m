classdef ProcessOfflineData < handle

    properties (Access = public)

    end
    properties (Access = private)
        mesh
        boundaryMeshJoined
        localGlobalConnecBd
        DirFun
        boundaryConditions
        bcApplier
        LHS
        RHS
        DDdofManager
        domainIndices

        fValuesTraining
        RigidBodyFun
        DeformationalFun

        fileNameData
        tolSameNode

    end


    properties (Access = private)
        nSubdomains
    end

    methods (Access = public)

        function obj = ProcessOfflineData()
            obj.init()
            obj.loadOfflineData();
%             data = load(obj.fileNameData);
%             obj.mesh            = data.mR;
%             obj.fValuesTraining = data.u;
            
            uFun   = obj.createDispFun();    
            RBfun  = obj.projectToRigidBody(uFun);
            DEFfun = obj.projectToDeformational(uFun,RBfun);

        end

    end

    methods (Access = private)

        function init(obj)
            obj.fileNameData = './Preconditioning/ROM/Training/PorousCell/OfflineData.mat';           
        end

        function loadOfflineData(obj)
            data = load(obj.fileNameData);
            obj.mesh            = data.mR;
            obj.fValuesTraining = data.u;
        end

        function uFun = createDispFun(obj)
            ntest  = size(obj.fValuesTraining,2);
            s.mesh  = obj.mesh;
            s.ndimf = 2;
            s.order = 'P1';
            for i = 1:ntest
                fValues   = obj.fValuesTraining(:,i);
                s.fValues = reshape(fValues,2,[])' ;
                uFun(i)   = LagrangianFunction(s);
            end
        end

        function RBfun = projectToRigidBody(obj,uFun)
            refPoint = (max(obj.mesh.coord)+min(obj.mesh.coord))/2;           
            ntest    = size(uFun,2);
            for i = 1:ntest
                RBfun(i) = uFun(i).project('RigidBody',refPoint);
            end
        end

        function DEFfun = projectToDeformational(obj,uFun,RBfun)         
            ntest    = size(uFun,2);
            
            for i = 1:ntest
                RBfun(i) = uFun(i).project('RigidBody',refPoint);
            end
        end

    end    

end
