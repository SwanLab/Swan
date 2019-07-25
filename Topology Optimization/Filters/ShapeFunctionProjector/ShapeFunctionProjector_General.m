classdef ShapeFunctionProjector_General < ShapeFunctionProjector
    
    properties (Access = private)
        unfittedMesh
        integrator 
        domainType
        interpolation
        quadrature        
    end
    
    methods (Access = public)
    
        function obj = ShapeFunctionProjector_General(cParams)
            obj.init(cParams)
            obj.domainType = cParams.domainType;
            obj.quadrature = cParams.quadrature;
            obj.createInterpolation();
            obj.createUnfittedMesh();
            obj.createIntegrator();                           
        end
    
        function xP = project(obj,x)
            nodalF = ones(size(x));
            obj.unfittedMesh.compute(x);
            %xP = obj.integrator.integrateUnfittedMesh(nodalF,obj.unfittedMesh);            
            xP = obj.integrator.integrate(nodalF);                        
        end        
        
    end
    
    methods (Access = private)
        
        function createInterpolation(obj)
            obj.interpolation = Interpolation.create(obj.mesh,'LINEAR');
            obj.interpolation.computeShapeDeriv(obj.quadrature.posgp)
        end      
        
        function createUnfittedMesh(obj)
            s.unfittedType = obj.domainType;
            s.meshBackground = obj.mesh;
            s.interpolationBackground = obj.interpolation;
            cParams = SettingsMeshUnfitted(s);
            obj.unfittedMesh = UnfittedMesh(cParams);            
        end
        
        function createIntegrator(obj)
            cParams.mesh = obj.unfittedMesh;
            cParams.type = obj.unfittedMesh.unfittedType;
            obj.integrator = Integrator.create(cParams);            
        end        
        
    end
    
end