classdef testComputingFemWithVademecumData < testShowingError ...
    & testLoadStoredVariable & testStoredComputedChecker
    
    properties (Access = protected)
        tol = 1e-6;
        testName = 'testComputingFemWithVademecumData';
        variablesToStore = {'u'};
    end
    
    properties (Access = private)
        fileName
        designVar
        Ctensor
        density
        fem
        densityPostProcess
    end
    
    methods (Access = public)
        
        function obj = testComputingFemWithVademecumData()
            obj.init();
            obj.createFEM();
            obj.createDesignVariableFromRandMxMy();
            obj.computeConstitutiveAndDensityFromVademecum();
            obj.computeFEM();
            obj.createDensityPostprocess();
            obj.printDensity();
            obj.selectComputedVar();
            
        end
    end
    
    methods (Access = protected)
        
        function selectComputedVar(obj)
            obj.computedVar{1} = obj.fem.variables.d_u;
        end

    end
    
    methods (Access = private)
        
        function init(obj)
            obj.fileName =  'VademecumSmoothCorner';
        end
        
        function createFEM(obj)
            obj.fem = FEM.create('test2d_quad');
            obj.fem.preProcess;
        end
        
        function createDesignVariableFromRandMxMy(obj)
            nelem = obj.fem.element.nelem;
            a = 0.01;
            b = 0.99;
           obj.designVar = (b-a).*[1:nelem;1:nelem]'/(nelem+1) + a;
        end
        
        function computeConstitutiveAndDensityFromVademecum(obj)
            s.fileName = obj.fileName;
            v = VademecumVariablesLoader(s);
            v.load();
            obj.Ctensor = v.Ctensor.compute(obj.designVar);
            obj.density = v.density.compute(obj.designVar);
        end
        
        function computeFEM(obj)
            obj.fem.setC(obj.Ctensor);
            obj.physicalProblem.computeStiffnessMatrix();
            obj.fem.computeVariables;
            obj.fem.print('ComputingFemVademecum');
        end
        
        function printDensity(obj)
            iter = 0;
            quad = Quadrature.set('QUAD');
            quad.computeQuadrature('CONSTANT');
            s.fields = obj.density;
            s.quad   = quad;
            obj.densityPostProcess.print(iter,s)
        end
        
        function createDensityPostprocess(obj)
            s.mesh    = obj.fem.mesh;
            s.outName = 'ComputingFemVademecum';
            ps = PostProcessDataBaseCreator(s);
            dB = ps.getValue();
            postCase = 'DensityGauss';
            obj.densityPostProcess = Postprocess(postCase,dB);
        end
        
    end
    
end