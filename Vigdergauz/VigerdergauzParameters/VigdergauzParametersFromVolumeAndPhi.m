classdef VigdergauzParametersFromVolumeAndPhi < handle
    
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
        
        function obj = VigdergauzParametersFromVolumeAndPhi(cParams)
            obj.init(cParams);
            obj.checkValidityOfParameters();
        end
        
        function compute(obj)
            r = obj.computeOptimalR();
            obj.computeAxAy(r)  
            obj.computeVigergauzParameters();
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.volume = cParams.volumeMicro;
            obj.phi    = atan(cParams.superEllipseRatio);
            obj.cx     = cParams.cx;
            obj.cy     = cParams.cy;            
        end
        
        function computeAxAy(obj,r)
            s.r = r;
            s.volume = obj.volume;
            s.cx = obj.cx;
            s.cy = obj.cy;
            axay = AxAyComputerFromVolumeAndR(s); 
            [obj.ax,obj.ay] = axay.compute();            
        end
        
        function r = computeOptimalR(obj)
            s.x0 = 1;
            s.functionToSolve = @(r) obj.equationForR(r);
            solver = ImplicitEquationSolver(s);
            r = solver.solve();
        end
               
        function f = equationForR(obj,r)
            obj.computeAxAy(r);
            obj.computeVigergauzParameters();
            mx = obj.parameters.mx;
            my = obj.parameters.my;
            f = tan(obj.phi) - mx/my;
        end
        
        function computeVigergauzParameters(obj)
            s.ax = obj.ax;
            s.ay = obj.ay;
            s.cx = obj.cx;
            s.cy = obj.cy;
            s.type = 'AxAndAy';
            obj.parameters = VigdergauzParameters.create(s);
        end 
        
        function checkValidityOfParameters(obj)
            if ~obj.isMicroStructureValid()
                error('Not possible axisRatio with this volume')
            end            
        end
        
        function itIs = isMicroStructureValid(obj)
            mMax = 0.99;
            phiMin = atan((1 - obj.volume)/(mMax^2));
            phiMax = atan((mMax^2)/(1 - obj.volume));
            itIs = obj.phi <= phiMax && obj.phi >= phiMin;
        end        
        
    end
    
end