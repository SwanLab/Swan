classdef Integrator_Composite < Integrator
    
    properties (GetAccess = public, SetAccess = private)
        integrators
        nInt
    end
    
    methods (Access = public)
        
        function obj = Integrator_Composite(cParams)
            obj.init(cParams);
            obj.createIntegrators(cParams);
        end
        
        function A = computeLHS(obj)
            npnod = obj.mesh.innerMeshOLD.npnod;
            A = sparse(npnod,npnod);
            for iInt = 1:obj.nInt
                A = A + obj.integrators{iInt}.computeLHS();
            end
        end
        
        function f = integrate(obj,nodalFunc)
            f = cell(1,obj.nInt);
            for iInt = 1:obj.nInt
                f{iInt} = obj.integrators{iInt}.integrate(nodalFunc);
            end
        end
        
    end
    
    methods (Access = private)
        
        function createNint(obj,cParams)
            obj.nInt = numel(cParams.compositeParams);
        end
        
        function createIntegrators(obj,cParams)
%             obj.createIntegratorsNew(cParams);
            obj.createIntegratorsOld();
        end
        
        function createIntegratorsNew(obj,cParams)
                        obj.createNint(cParams);
            params = cParams.compositeParams;
            for iInt = 1:obj.nInt
                s = params{iInt};
                integrator = Integrator.create(s);
                obj.integrators{end+1} = integrator;
            end
        end
        
        function createIntegratorsOld(obj)
            meshC = obj.mesh;
            activeMeshes = meshC.getActiveMeshes();
            for iMesh = 1:meshC.nActiveMeshes
                thisMesh = activeMeshes{iMesh};
                cParams.mesh = activeMeshes{iMesh};
                cParams.type = activeMeshes{iMesh}.unfittedType;
                switch thisMesh.unfittedType
                    case 'SIMPLE'
                        cParams.globalConnec = meshC.globalConnectivities{iMesh};
                        cParams.backgroundMesh = thisMesh;
                        cParams.innerToBackground = [];
                        cParams.npnod = obj.mesh.innerMeshOLD.npnod;
                        obj.integrators{iMesh} = Integrator.create(cParams);
                    case {'INTERIOR','BOUNDARY'}
                        obj.integrators{iMesh} = Integrator.create(cParams);
                end
            end
        end
        
    end
    
end

