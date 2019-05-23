classdef Pedersen < handle
    
    properties (GetAccess = public, SetAccess = private)
        optimalAngle
        compliance
        RotationMatrix        
    end
    
    properties (Access = private)
        angle
        strain
        Ctensor
    end
    
    methods (Access = public)
        
        function obj = Pedersen()
            obj.init();
        end
        
        function compute(obj,strain,Ctensor)
            obj.strain = strain;
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
            R = [ 0.5*(1+c) 0.5*(1-c) s;
                0.5*(1-c) 0.5*(1+c) -s;
                -s/2       s/2      c];
            obj.RotationMatrix = R;
        end
        
        function computeCompliance(obj)
            e = obj.strain;
            C = obj.Ctensor;
            R = obj.RotationMatrix;
            c = e'*R*C*R'*e;
            %Re = inv(R');
            %c = e'*Re'*C*Re*e;            
            obj.compliance = matlabFunction(c);
        end
        
        function obtainOptimalAngle(obj)
           x = fminbnd(obj.compliance,0,2*pi);
           obj.optimalAngle = x;
        end
        
        
    end
end
