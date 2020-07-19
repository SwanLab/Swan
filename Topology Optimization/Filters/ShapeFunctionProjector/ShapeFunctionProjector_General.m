classdef ShapeFunctionProjector_General < ShapeFunctionProjector
    
    properties (Access = private)
        unfittedMesh
        domainType
    end
    
    methods (Access = public)
        
        function obj = ShapeFunctionProjector_General(cParams)
            obj.init(cParams)
            obj.domainType = cParams.domainType;
            obj.createUnfittedMesh();
        end
        
        function fInt = project(obj,ls)
            if all(ls>0)
                fInt = zeros(size(ls));
            else
                fNodes = ones(size(ls));
                obj.unfittedMesh.compute(ls);    
%                obj.unfittedMesh.plotBoundary(); 
%                view([1 1 1]);
%                drawnow
                s.mesh = obj.unfittedMesh;
                s.type = 'Unfitted';
                integrator = Integrator.create(s);            
                fInt = integrator.integrateInDomain(fNodes);                
            end
            
        end
        
    end
    
    methods (Access = private)
        
        
        function createUnfittedMesh(obj)
            s.backgroundMesh = obj.mesh.innerMeshOLD;
            s.boundaryMesh   = obj.mesh.boxFaceMeshes;
            cParams = SettingsMeshUnfitted(s);
            obj.unfittedMesh = UnfittedMesh(cParams);
        end
        
    end
    
end