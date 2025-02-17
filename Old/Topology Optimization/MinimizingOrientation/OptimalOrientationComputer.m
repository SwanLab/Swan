classdef OptimalOrientationComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        optimalAngle
        compliance
        Rsigma
    end
    
    properties (Access = private)
        angle
        stress
        Ctensor
    end
    
    methods (Access = public)
        
        function obj = OptimalOrientationComputer()
            obj.init();
        end
             
        function compute(obj,stress,Ctensor)
            obj.stress = stress;
            obj.Ctensor = Ctensor;
            obj.computeCompliance();
            obj.obtainOptimalAngle();
        end
       
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.createAngle();
            obj.createRotationMatrix();
        end
        
        function createAngle(obj)
            obj.angle = sym('a','real');
        end
        
        function createRotationMatrix(obj)
            a = obj.angle;
            c = cos(2*a);
            s = sin(2*a);
            Rs = [ 0.5*(1+c) 0.5*(1-c) s;
                 0.5*(1-c) 0.5*(1+c) -s;
                 -s/2       s/2      c];
            obj.Rsigma = Rs;
        end
        
        function computeCompliance(obj)
            s = obj.stress;
            C = obj.Ctensor;
            Rs = obj.Rsigma;
            c = s'*Rs'*inv(C)*Rs*s;
            obj.compliance = matlabFunction(c);
        end
        
        function obtainOptimalAngle(obj)
           x = fminbnd(obj.compliance,-pi,pi);
           obj.optimalAngle = x;
        end
        
    end
end
