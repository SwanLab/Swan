classdef Integrator_Composite < Integrator
    
    properties (GetAccess = public, SetAccess = private)
        integrators
        nInt
    end
    
    methods (Access = public)
        
        function obj = Integrator_Composite(mesh)
            obj.createIntegrators(mesh);
        end
        
    end
    
    methods (Access = protected)
        
        function f = computeIntegral(obj,nodalFunc)
            f = cell(1,obj.nInt);
            for iInt = 1:obj.nInt
                f{iInt} = obj.integrators{iInt}.computeIntegral(nodalFunc);
            end
        end
        
        function A = computeLHS(obj)
            npnod = obj.meshBackground.npnod;
            A = sparse(npnod,npnod);
            for iInt = 1:obj.nInt
                globalConnec = obj.meshUnfitted.globalConnectivities{iInt};
                A = A + obj.integrators{iInt}.computeLHS(globalConnec,npnod);
            end
        end
        
    end
    
    methods (Access = private)
        
        function createIntegrators(obj,meshComposite)
            activeMeshes = meshComposite.getActiveMeshes();
            for iMesh = 1:meshComposite.nActiveMeshes
                mesh = activeMeshes{iMesh};
                obj.integrators{iMesh} = Integrator.create(mesh);
            end
            obj.nInt = numel(obj.integrators);
        end
                
    end
    
end

