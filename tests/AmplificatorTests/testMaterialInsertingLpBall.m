classdef testMaterialInsertingLpBall < testShowingError
    
    properties (Access = protected)
        tol = 5e-2;
        testName = 'testMaterialInsertingLpBall';
    end
    
    properties (Access = private)  

    end
    
    
    methods (Access = public)
        
        function obj = testMaterialInsertingLpBall()
            matProp.rho_plus = 1;
            matProp.rho_minus = 0;
            matProp.E_plus = 1;
            matProp.E_minus = 1e-3;
            matProp.nu_plus = 1/3;
            matProp.nu_minus = 1/3;
            d.constitutiveProperties = matProp;
            d.interpolation = 'SIMPALL';
            d.dim = '2D';
            d.typeOfMaterial = 'ISOTROPIC';
            mI = Material_Interpolation.create(d);
            mI.computeMatProp([0,1])
        end
        
    end
    
    methods (Access = protected)
        
        function computeError(obj)
            obj.error = 0;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.pNorm = 8;
        end
        

        
    end
    
    
end