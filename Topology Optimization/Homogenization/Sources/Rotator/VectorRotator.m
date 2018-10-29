classdef VectorRotator < Rotator
    
    properties (Access = protected)
        rotatedTensor
        rotationMatrix
    end
    
    properties (Access = private)
        dim
        FirstTerm
        SecondTerm
        ThirdTerm
    end
    
    methods (Access = public)
        
        function obj = VectorRotator(angle,dir)
            obj.init(angle,dir)
            obj.generateRotator()
        end
        
    end
    
    methods (Access = protected)
        function computeRotation(obj,vector)
            R = obj.rotationMatrix;
            s = vector;
            obj.rotatedTensor = R*s;
        end
        
        function init(obj,angle,vect)
            obj.init@Rotator(angle,vect)
            obj.dim = 3;
        end
    end
    
    methods (Access = private)
        
        function generateRotator(obj)
            obj.computeTerms()
            obj.addTerms()
        end
        
        function computeTerms(obj)
            obj.computeFirstTerm()
            obj.computeSecondTerm()
            obj.computeThirdTerm()
        end
        
        function computeFirstTerm(obj)
            I = eye(obj.dim);
            alpha = obj.angle;
            obj.FirstTerm = cos(alpha)*I;
        end
        
        function computeSecondTerm(obj)
            u(:,1) = obj.dir;
            alpha = obj.angle;
            A = sym(zeros(obj.dim,obj.dim));
            A(1,2) = -u(3,1);
            A(1,3) = u(2,1);
            A(2,1) = u(3,1);
            A(2,3) = -u(1,1);
            A(3,1) = -u(2,1);
            A(3,2) = u(1,1);
            obj.SecondTerm = sin(alpha)*A;
        end
        
        function computeThirdTerm(obj)
            u(:,1) = obj.dir;
            alpha = obj.angle;
            obj.ThirdTerm = (1-cos(alpha))*(u*u');
        end
        
        function addTerms(obj)
            R1 = obj.FirstTerm;
            R2 = obj.SecondTerm;
            R3 = obj.ThirdTerm;
            obj.rotationMatrix = R1 + R2 + R3;
        end
        
    end

end

