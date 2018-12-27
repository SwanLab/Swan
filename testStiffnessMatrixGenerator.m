classdef testStiffnessMatrixGenerator < ...
        testShowingError & ...
        testLoadStoredVariable & ...
        testStoredComputedChecker
    
    properties (Access = protected)
        testName = 'test2dFourquad';
        variablesToStore = {'K','KwithVoidMaterial'};
        tol = 1e-6;
    end
    
    properties (Access = private)
        MatInterp
        density
        voidDens
        fem
        femVoid
    end
    
    methods (Access = public)
        
        function obj = testStiffnessMatrixGenerator()
            obj.createMaterialInterpolation()
            obj.createFemProblems()
            obj.createDensities()
            obj.computeFemProblems()
            obj.selectComputedVar()
        end
        
    end
    
    methods (Access = protected)
        
        function selectComputedVar(obj)
            obj.computedVar{1} = obj.fem.element.K(:);
            obj.computedVar{2} = obj.femVoid.element.K(:); 
        end
        
    end
    
    methods (Access = private)
        
        function createMaterialInterpolation(obj)
            matProp.rho_plus = 1;
            matProp.rho_minus = 0;
            matProp.E_plus = 1;
            matProp.E_minus = 1e-3;
            matProp.nu_plus = 1/3;
            matProp.nu_minus = 1/3;
            obj.MatInterp = Material_Interpolation.create(matProp,'ISOTROPIC','SIMPALL','2D');
        end
        
        function createFemProblems(obj)
           obj.fem = FEM.create(obj.testName);
           obj.femVoid = FEM.create(obj.testName);           
        end

        function createDensities(obj)
            obj.density = ones(obj.fem.geometry.interpolation.nelem,1);
            obj.voidDens = obj.density;
            obj.voidDens([1 4]) = 0.01;
        end
        
        function computeFemProblems(obj)
            obj.computeFemProblem(obj.density,obj.fem)
            obj.computeFemProblem(obj.voidDens,obj.femVoid)            
        end
        
        function computeFemProblem(obj,dens,fem)
            matProp = obj.MatInterp.computeMatProp(dens);
            fem.preProcess();
            fem.setMatProps(matProp);
            fem.computeVariables();
        end
        
    end
    
    
end

