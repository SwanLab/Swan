classdef Integrator_Composite < Integrator
    
    properties (GetAccess = public, SetAccess = private)
        integrators
        nInt
    end
    
    properties (Access = private)
        integratorInterior
        integratorsBoxFaces
    end
    
    methods (Access = public)
        
        function obj = Integrator_Composite(mesh)
            obj.createIntegrators(mesh);
            %             obj.createInteriorIntegrator(meshComposite);
            %             obj.createBoxFacesIntegrators(meshComposite);
        end
        
    end
    
    methods (Access = protected)
        
        function f = computeIntegral(obj,nodalFunc)
            f = cell(1,obj.nInt);
            for iInt = 1:obj.nInt
                f{iInt} = obj.integrators{iInt}.computeIntegral(nodalFunc);
            end
            %             A.interiorIntegral = obj.computeInteriorIntegral(F1);
            %             A.boxFacesIntegrals = obj.computeBoxFacesIntegrals(F1);
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
        
        function createInteriorIntegrator(obj,meshComposite)
            obj.integratorInterior = Integrator.create(meshComposite.meshInterior);
        end
        
        function createBoxFacesIntegrators(obj,meshComposite)
            for iactive = 1:meshComposite.nActiveBoxFaces
                iface = meshComposite.activeBoxFaceMeshesList(iactive);
                obj.integratorsBoxFaces{iface} = Integrator.create(meshComposite.boxFaceMeshes{iface});
            end
        end
        
        function A = computeInteriorIntegral(obj,F1)
            A = obj.integratorInterior.computeIntegral(F1);
        end
        
        function A = computeBoxFacesIntegrals(obj,F1)
            A = cell(1,obj.meshUnfitted.nActiveBoxFaces);
            for iactive = 1:obj.meshUnfitted.nActiveBoxFaces
                iface = obj.meshUnfitted.activeBoxFaceMeshesList(iactive);
                A{iface} = obj.integratorsBoxFaces{iface}.computeIntegral(F1);
            end
        end
        
    end
    
end

