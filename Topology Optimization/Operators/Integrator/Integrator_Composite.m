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
            npnod = obj.mesh.innerMesh.npnod;
            A = sparse(npnod,npnod);
            for iInt = 1:obj.nInt
                globalConnec = obj.mesh.globalConnectivities{iInt};
                A = A + obj.integrators{iInt}.computeLHS(globalConnec,npnod);
            end
        end 
        
        function f = computeIntegral(obj,nodalFunc)
            f = cell(1,obj.nInt);
            for iInt = 1:obj.nInt
                f{iInt} = obj.integrators{iInt}.computeIntegral(nodalFunc);
            end
        end        
        
    end
    
    methods (Access = protected)
        

        

        
    end
    
    methods (Access = private)
        
        function createIntegrators(obj)
            meshC = obj.mesh;
            activeMeshes = meshC.getActiveMeshes();
            for iMesh = 1:meshC.nActiveMeshes
                cParams.mesh = activeMeshes{iMesh};                
                obj.integrators{iMesh} = Integrator.create(cParams);
            end
            obj.nInt = numel(obj.integrators);
        end
                
    end
    
end

