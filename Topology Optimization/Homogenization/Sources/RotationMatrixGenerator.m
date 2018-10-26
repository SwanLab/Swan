classdef RotationMatrixGenerator < handle
    
    
    properties (Access = public)
        Matrix
    end
    
    properties (Access = private)
        Alpha
        NormalVector
        Dim
        FirstTerm
        SecondTerm
        ThirdTerm
    end
    
    methods (Access = public)
        
        function obj = RotationMatrixGenerator(Angle,NormalVector)            
            obj.init(Angle,NormalVector);
            obj.computeMatrix()
        end
        
        function r = getValue(obj)
            r = double(obj.Matrix);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,Angle,Vect)
            obj.Alpha = Angle;
            obj.NormalVector = Vect;    
            obj.Dim = 3;
        end
        
        function computeMatrix(obj)
            obj.computeTerms()
            obj.addTerms()
        end
        
        function computeTerms(obj)
            obj.computeFirstTerm()
            obj.computeSecondTerm()
            obj.computeThirdTerm()            
        end
        
        function computeFirstTerm(obj)
            I = eye(obj.Dim);
            alpha = obj.Alpha;            
            obj.FirstTerm = cos(alpha)*I;            
        end
        
        function computeSecondTerm(obj)
            u(:,1) = obj.NormalVector;
            alpha = obj.Alpha; 
            A = sym(zeros(obj.Dim,obj.Dim));
            A(1,2) = -u(3,1);
            A(1,3) = u(2,1);
            A(2,1) = u(3,1);
            A(2,3) = -u(1,1);
            A(3,1) = -u(2,1);
            A(3,2) = u(1,1);            
            obj.SecondTerm = sin(alpha)*A;
        end
        
        function computeThirdTerm(obj)
            u(:,1) = obj.NormalVector;
            alpha = obj.Alpha;
            obj.ThirdTerm = (1-cos(alpha))*(u*u');
        end
        
        function addTerms(obj)
            R1 = obj.FirstTerm;
            R2 = obj.SecondTerm;
            R3 = obj.ThirdTerm;
            obj.Matrix = R1 + R2 + R3;
        end
        
    end
    
end

