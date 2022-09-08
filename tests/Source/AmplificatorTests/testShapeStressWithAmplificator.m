classdef testShapeStressWithAmplificator < testShowingError
    
    properties (Access = protected)
        tol = 5e-2;
        testName = 'testShapeStressWithAmplificator';
    end
    
    properties (Access = private)
        PcompHomog
        PcompShape
        m1 = 0.2;
        m2 = 0.8;
        printing = true;
        iter = 0;    
        homog
        strain
    end
    
    
    methods (Access = public)
        
        function obj = testShapeStressWithAmplificator()
            obj.createSmoothRectInclusionHomogenizer()
            obj.obtainPcomponentWithHomogenizer(1,1);
            obj.obtainPcomponentWithShapeFunction();
        end
        
    end
    
    methods (Access = protected)
        
        function computeError(obj)
            obj.error = abs(obj.PcompHomog - obj.PcompShape)/abs(obj.PcompShape);
        end
        
    end
    
    methods (Access = private)
        
        function createSmoothRectInclusionHomogenizer(obj)
            f    = strcat(obj.testName);
            p    = obj.printing;
            m1v  = obj.m1;
            m2v  = obj.m2;
            i    = obj.iter;
            obj.homog = NumericalSmoothRectangleHomogenizer(f,p,m1v,m2v,i);
        end
        
        function obtainPcomponentWithHomogenizer(obj,i,j)
            Ch            = obj.homog.getCh();
            Ptensor       = obj.homog.getAmplificatorTensor();
            PTensorValues = Ptensor.getValue();
            obj.PcompHomog = PTensorValues(i,j);
        end
        
        function computeStressShape(obj)
            obj.computeStrainWithCanonicalStress();
            settings = obj.homog.getSettings();
            ls = obj.homog.getLevelSet();
            sF = ShFunc_StressNorm(settings);
%             sF.filter.preProcess();
            sF.setVstrain(obj.strain);
            sF.computeCostAndGradient(ls);
            obj.stressShape = sF;
        end
        
        function obtainPcomponentWithShapeFunction(obj)
            obj.PcompShape = obj.stressShape.getValue();
        end        
        
        function storeMicroProblem(obj)
            obj.microProblem = obj.stressShape.getPhysicalProblems();
        end
        
        function computeStrainWithCanonicalStress(obj)
            stress = [1,0,0]';
            Ch = obj.homog.getCh();           
            obj.strain = Ch\stress;
        end
        
    end
    
    
end