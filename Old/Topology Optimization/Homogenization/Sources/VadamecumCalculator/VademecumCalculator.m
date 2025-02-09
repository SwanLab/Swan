classdef VademecumCalculator < handle
    
    properties (Access = public)
        cellVariables        
    end
    
    properties (Access = private)
        optimalExponent
        qOptimal        
    end
    
    properties (Access = private)            
        rho
        xi
        phi 
        nCase
    end
    
    methods (Access = public)
        
        function obj = VademecumCalculator(cParams)
            obj.init(cParams)
            obj.createOptimalExponentComputer();            
        end
        
        function compute(obj)
            obj.computeOptimalExponent();
            obj.obtainCellVariables(obj.qOptimal);
        end
        
        function obtainCellVariables(obj,q)
            cVariables = obj.optimalExponent.obtainCellVariables(q);  
            cVariables.rho = obj.rho;
            cVariables.xi = obj.xi;
            cVariables.phi = obj.phi;
            cVariables.q = q;
            obj.cellVariables = cVariables;
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.xi  = cParams.xi;
            obj.rho = cParams.rho;
            obj.phi = cParams.phi;       
            obj.nCase = cParams.nCase;
        end
        
        function createOptimalExponentComputer(obj)
            s.fileName = ['OptimaSuperEllipseCase',obj.nCase];
            s.rho   = obj.rho;
            s.txi   = obj.xi;
            s.phi   = obj.phi;
            s.pNorm = 'max';
            s.hMesh = [];            
            optimalExp = OneOptimalExponentComputerAndFunctionVariation(s);           
            obj.optimalExponent = optimalExp;            
        end        
        
        function q = computeOptimalExponent(obj)
             obj.optimalExponent.computeOptimalExponent();  
             qMin = obj.optimalExponent.qOptIter();
             fMin = obj.optimalExponent.fOptIter();            
             [~,ind] = min(fMin);
             q = qMin(ind);
             obj.qOptimal = q;
        end
        

        
    end
    
    
    
end