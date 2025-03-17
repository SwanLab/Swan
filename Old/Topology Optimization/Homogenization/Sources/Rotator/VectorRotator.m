classdef VectorRotator < Rotator
    
    properties (Access = private)
        dim
        FirstTerm
        SecondTerm
        ThirdTerm
    end
    
    methods (Access = public)
        
        function obj = VectorRotator(angle,dir)
            obj.compute(angle,dir)
        end
        
    end
    
    methods (Access = protected)
        function computeRotation(obj,vector)
            obj.createRotatedTensor(vector);
            R = obj.rotationMatrix;
            s = vector.getValue;
            sR = R*s;
            obj.rotatedTensor.setValue(sR);
        end
        
        function init(obj,angle,vect)
            obj.init@Rotator(angle,vect)
            obj.dim = 3;
        end
        
        function generateRotator(obj)
            obj.computeTerms()
            obj.addTerms()
        end

        
    end
    
    methods (Access = private)
        
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
            u(:,1) = obj.dir.getValue();
            alpha = obj.angle;
            A = obj.createSecondTermMatrix();
            A(1,2) = -u(3,1);
            A(1,3) = u(2,1);
            A(2,1) = u(3,1);
            A(2,3) = -u(1,1);
            A(3,1) = -u(2,1);
            A(3,2) = u(1,1);
            obj.SecondTerm = sin(alpha)*A;
        end
        
        function computeThirdTerm(obj)
            u(:,1) = obj.dir.getValue();
            alpha = obj.angle;
            obj.ThirdTerm = (1-cos(alpha))*(u*u');
        end
        
        function addTerms(obj)
            R1 = obj.FirstTerm;
            R2 = obj.SecondTerm;
            R3 = obj.ThirdTerm;
            obj.rotationMatrix = R1 + R2 + R3;
        end
        
        function A = createSecondTermMatrix(obj)
            if obj.isSymbolic()
              A = sym(zeros(obj.dim,obj.dim));
            else
              A = zeros(obj.dim,obj.dim);
            end
        end
        
        function itIs = isSymbolic(obj)
            isDirSym = isa(obj.dir.getValue(),'sym');
            isAngleSym = isa(obj.angle,'sym');
            itIs = isDirSym || isAngleSym;
        end
        
    end

end

