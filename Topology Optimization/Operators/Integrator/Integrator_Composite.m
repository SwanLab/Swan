classdef Integrator_Composite < Integrator
    
    properties (Access = private)
        integratorInterior
        integratorsBoxFaces
    end
    
    methods (Access = public)
        
        function obj = Integrator_Composite(meshComposite)
            obj.createInteriorIntegrator(meshComposite);
            obj.createBoxFacesIntegrators(meshComposite);
        end
        
    end
    
    methods (Access = protected)
        
        function A = computeIntegral(obj,F1)
            A.interiorIntegral = obj.computeInteriorIntegral(F1);
            A.boxFacesIntegrals = obj.computeBoxFacesIntegrals(F1);
        end
        
    end
    
    methods (Access = private)
        
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

