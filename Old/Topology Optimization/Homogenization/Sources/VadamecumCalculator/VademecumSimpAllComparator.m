classdef VademecumSimpAllComparator < handle
    
    properties (Access = public)
        density
        CtensorVademecum
        CtensorSIMPALL        
    end
    
    properties (Access = private)
        fileName
        designVar
    end
    
    methods (Access = public)
        
        function obj = VademecumSimpAllComparator(cParams)
            obj.init(cParams);
        end
        
        function calculate(obj)
            obj.createDesignVariableFromRandMxMy();
            obj.computeConstitutiveAndDensityFromVademecum();            
            obj.computeConstitutiveFromDensity();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = cParams.fileName;
        end
        
        function computeConstitutiveAndDensityFromVademecum(obj)
            s.fileName = obj.fileName;
            v = VademecumVariablesLoader(s);
            v.load();
            obj.CtensorVademecum = v.Ctensor.compute(obj.designVar);
            obj.density = v.density.compute(obj.designVar);
        end
        
        function createDesignVariableFromRandMxMy(obj)
            a = 0.01;
            b = 0.99;
            obj.designVar = (b-a).*rand(1000,2) + a;
        end
        
        
        function computeConstitutiveFromDensity(obj)
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
            material = mI.computeMatProp(obj.density);
            me = Isotropic2dElasticMaterial(size(obj.density,1));
            me.setProps(material);
            obj.CtensorSIMPALL = me.C;
        end
        
        
    end
    
    
end

