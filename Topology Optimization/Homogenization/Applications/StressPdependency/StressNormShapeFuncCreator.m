classdef StressNormShapeFuncCreator < handle
    
    properties (Access = private)
        outFile
        stressShape
        stressShapeDB
    end
    
    methods (Access = public)
        
        function obj = StressNormShapeFuncCreator(d)
            obj.init(d)
            obj.createStressShapeDataBase();
            obj.createStressShape();
        end
        
        function s = getStressNormShape(obj)
            s = obj.stressShape;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,d)
            obj.outFile       = d.outFile;
        end
        
        function createStressShape(obj)
            dB = obj.stressShapeDB();
            sF = ShFunc_StressNorm(dB);
%             sF.filter.preProcess();
            obj.stressShape = sF;
        end
        
        function createStressShapeDataBase(obj)
            dB.filename    = obj.outFile;            
            dB.TOL         = obj.createMaterialProperties();
            dB.material    = 'ISOTROPIC';
            dB.method      = 'SIMPALL';
            dB.pdim        = '2D';
            dB.stressHomog = [1,0,0]';
            dB.filter      = 'P1';
            dB.optimizer   = 'SLERP';
            dB.pdim        = '2D';
            dB.scale       = 'MICRO';
            obj.stressShapeDB = dB;
        end
    end
    
    methods (Access = private, Static)
        
        function matProp = createMaterialProperties()
            matProp.rho_plus  =  1;
            matProp.rho_minus =  0;
            matProp.E_plus    =  1;
            matProp.E_minus   =  1.0000e-03;
            matProp.nu_plus   =  0.3333;
            matProp.nu_minus  =  0.3333;
        end
        
    end
    
end