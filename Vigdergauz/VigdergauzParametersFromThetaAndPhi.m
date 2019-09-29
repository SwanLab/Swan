classdef VigdergauzParametersFromThetaAndPhi < handle
    
    properties (Access = public)
        parameters
    end
    
    properties (Access = private)
        cx
        cy
        phi
        volume
        ax
        ay        
    end
    
    methods (Access = public)
        
        function obj = VigdergauzParametersFromThetaAndPhi(cParams)
            obj.init(cParams);
            obj.checkValidityOfParameters();
            obj.compute();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.volume = cParams.volume;
            obj.phi    = cParams.phi;
            obj.cx     = cParams.cx;
            obj.cy     = cParams.cy;            
        end
        
        function compute(obj)
            r = obj.computeOptimalR();
            axay = AxAyComputerFromVolumeAndR(obj.volume,r,obj.cx);
            [obj.ax,obj.ay] = axay.compute();  
            obj.computeVigergauzParameters();
        end                
        
        function r = computeOptimalR(obj)
            s.x0 = 1;
            s.functionToSolve = @(r) obj.equationForR(r);
            solver = ImplicitEquationSolver(s);
            r = solver.solve();
        end
               
        function f = equationForR(obj,r)
            axay = AxAyComputerFromVolumeAndR(obj.volume,r,obj.cx);
            [obj.ax,obj.ay] = axay.compute();
            obj.computeVigergauzParameters();
            rx = obj.parameters.rx;
            ry = obj.parameters.ry;
            f = tan(obj.phi) - rx/ry;
        end
        
        function computeVigergauzParameters(obj)
            s.ax = obj.ax;
            s.ay = obj.ay;
            s.cx = obj.cx;
            s.cy = obj.cy;
            obj.parameters = VigdergauzParametersComputerFromAxAy(s);
        end 
        
        function checkValidityOfParameters(obj)
            if ~obj.isMicroStructureValid()
                error('Not possible axisRatio with this volume')
            end            
        end
        
        function itIs = isMicroStructureValid(obj)
            rMax = 0.99;
            phiMin = atan((obj.volume)/(rMax^2));
            phiMax = atan((rMax^2)/(obj.volume));
            itIs = obj.phi <= phiMax && obj.phi >= phiMin;
        end        
        
    end
    
end