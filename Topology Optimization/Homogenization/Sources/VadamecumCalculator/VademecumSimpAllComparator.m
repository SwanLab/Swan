classdef VademecumSimpAllComparator < handle
    
    properties (Access = public)
        density
        CtensorVademecum
        CtensorSIMPALL        
    end
    
    properties (Access = private)
        fileName
        vadVariables
        interpolator
        designVar
    end
    
    methods (Access = public)
        
        function obj = VademecumSimpAllComparator(cParams)
            obj.init(cParams);
        end
        
        function calculate(obj)
            obj.createDesignVariableFromRandMxMy();
            obj.loadVademecumVariables();
            obj.createInterpolator();
            obj.computeConstitutiveFromVademecum();
            obj.computeDensityFromVademecum();
            obj.computeConstitutiveFromDensity();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = cParams.fileName;
        end
        
        function createDesignVariableFromRandMxMy(obj)
            a = 0.01;
            b = 0.99;
            obj.designVar = (b-a).*rand(1000,2) + a;
        end
        
        function loadVademecumVariables(obj)
            matFile   = [obj.fileName,'.mat'];
            file2load = fullfile('Output',obj.fileName,matFile);
            v = load(file2load);
            obj.vadVariables = v.d;
        end
        
        function createInterpolator(obj)
            sM.x = obj.vadVariables.domVariables.mxV;
            sM.y = obj.vadVariables.domVariables.myV;
            sI.mesh = StructuredMesh(sM);
            obj.interpolator = Interpolator(sI);
        end
        
        function computeConstitutiveFromVademecum(obj)
            s.vadVariables = obj.vadVariables;
            s.interpolator = obj.interpolator;
            ct = ConstitutiveTensorFromVademecum(s);
            obj.CtensorVademecum = ct.computeCtensor(obj.designVar);
        end
        
        function computeDensityFromVademecum(obj)
            s.vadVariables = obj.vadVariables;
            s.interpolator = obj.interpolator;
            dt = DensityFromVademecum(s);
            obj.density = dt.computeDensity(obj.designVar);
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
            me = Material_Elastic_ISO_2D(size(obj.density,1));
            me.setProps(material);
            obj.CtensorSIMPALL = me.C;
        end
        
        
    end
    
    
end

