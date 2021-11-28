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
                s.mesh = obj.unfittedMesh;
                s.type = 'Unfitted';
                s.dim  = obj.computeDim(s.mesh.backgroundMesh);
                s.npnod = obj.unfittedMesh.backgroundMesh.npnod;
                integrator = Integrator.create(s);
                fInt = integrator.integrateInDomain(fNodes);
            end
            
        end
        
    end
    
    methods (Access = private)
        
        function createUnfittedMesh(obj)
            s.backgroundMesh = obj.mesh;
            cParams = SettingsMeshUnfitted(s);
            obj.unfittedMesh = UnfittedMesh(cParams);
        end
        
        function dim = computeDim(obj, mesh)
            s.ngaus = [];
            s.mesh  = obj.mesh;
            s.pdim  = obj.createPdim(mesh);
            dim    = DimensionVariables(s);
            dim.compute();
        end

        function pdim = createPdim(obj, mesh)
            switch mesh.ndim
                case 2
                    pdim = '2D';
                case 3
                    pdim = '3D';
            end
        end

    end
    
end