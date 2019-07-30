classdef Integrator_Composite < Integrator
    
    properties (GetAccess = public, SetAccess = private)
        integrators
        nInt
    end
    
    methods (Access = public)
        
        function obj = Integrator_Composite(cParams)
            obj.init(cParams);
            obj.createIntegrators();
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
        
        function createIntegrators(obj)
           
%             cParamsInnerCut.mesh = obj.mesh.innerCutMesh;
%             cParamsInnerCut.type = 'CutMesh';
%             innerCutIntegrator = Integrator.create(cParamsInnerCut);
%             obj.integrators{end+1} = innerCutIntegrator;
%             
%             cParamsInner.mesh = obj.mesh.innerMesh;
%             cParamsInner.type = 'SIMPLE';
%             cParamsInner.globalConnec = obj.mesh.globalConnec;
%             cParamsInner.npnod = obj.mesh.innerMesh.npnod;
%             innerIntegrator = Integrator.create(cParamsInner);
%             obj.integrators{end+1} = innerIntegrator;

            
            meshC = obj.mesh;
            activeMeshes = meshC.getActiveMeshes();
            for iMesh = 1:meshC.nActiveMeshes
                thisMesh = activeMeshes{iMesh};
                cParams.mesh = activeMeshes{iMesh};
                cParams.type = activeMeshes{iMesh}.unfittedType;
                switch thisMesh.unfittedType
                    case 'SIMPLE'
                        cParams.globalConnec = meshC.globalConnectivities{iMesh};
                        cParams.npnod = obj.mesh.innerMeshOLD.npnod;
                        obj.integrators{iMesh} = Integrator.create(cParams);
                    case {'INTERIOR','BOUNDARY'}
                        obj.integrators{iMesh} = Integrator.create(cParams);
                end
            end
        end
        
    end
    
    methods
        
        function n = get.nInt(obj)
            n = numel(obj.integrators);
        end
        
    end
    
end

