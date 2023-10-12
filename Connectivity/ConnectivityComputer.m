classdef ConnectivityComputer < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
       mesh
       levelSet
       density
       conductivity
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = ConnectivityComputer()
            obj.init();
            obj.createMesh();
            obj.createLevelSet();
            obj.filterCharacteristicFunction();
            obj.createMaterialInterpolator();
            obj.createStiffnessMatrix();
            obj.density.plot()
            figure
            obj.levelSet.getUnfittedMesh().plot()
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            
        end
        
        function createMesh(obj)
            x1 = linspace(-1,1,20);
            x2 = linspace(0,1,20);
            [xv,yv] = meshgrid(x1,x2);
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            m = Mesh(s);            
            obj.mesh = m;
        end

        function createLevelSet(obj)
            s.ndim       = 2;
            s.fracRadius = 0.4;
            s.coord      = obj.mesh.coord;
            sD.type = 'LevelSet';
            sD.mesh = obj.mesh;
            sD.creatorSettings = s;
            sD.initialCase = 'circleInclusion';
            obj.levelSet   = DesignVariable.create(sD);
        end
 
        function filterCharacteristicFunction(obj)
            s.mesh = obj.mesh;
            s.quadratureOrder = [];
            s.femSettings.scale = 'MACRO';
            s.designVariable = obj.levelSet;
            f = Filter_PDE_LevelSet(s);
            dens = f.getP0fromP1([]);
            s.fValues = dens;
            s.mesh    = obj.mesh;
            obj.density = P0Function(s);
        end

        function createMaterialInterpolator(obj)
            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPThermal';
            s.alpha0         = 1e-3;
            s.alpha1         = 1;
            m = MaterialInterpolation.create(s);
            mp = m.computeMatProp(obj.density.fValues);
            s.fHandle = mp.alpha;
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            a = AnalyticalFunction(s);
            obj.conductivity = a;
        end

        function createStiffnessMatrix(obj)
            s.test  = P1Function.create(obj.mesh,1); 
            s.trial = P1Function.create(obj.mesh,1); 
            s.mesh  = obj.mesh;
            s.quadratureOrder = 'QUADRATIC';
            s.function        = obj.conductivity;
            s.type            = 'StiffnessMatrixWithFunction';
            LHS = LHSintegrator.create(s);
        end
        
    end
    
end