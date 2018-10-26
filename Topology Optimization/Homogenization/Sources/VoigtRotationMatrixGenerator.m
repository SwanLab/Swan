classdef VoigtRotationMatrixGenerator < handle
    
    properties (Access = public)
        VoigtMatrix
        VoigtMatrixPlaneStress
    end

    
    properties (Access = private)
        Angle
        Direction
    end
    
    
    methods (Access = public)
        
        function obj = VoigtRotationMatrixGenerator(Angle,Direction)
            obj.init(Angle,Direction)
            obj.createVoigtMatrix()
            obj.createVoigtMatrixPlaneStress()
        end
        
        function r = getValue(obj)
            r = double(obj.VoigtMatrix);
        end
        
        function r = getPsValue(obj)
            r = double(obj.VoigtMatrixPlaneStress);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,Angle,Direction)
            obj.Angle = Angle;
            obj.Direction = Direction;
        end
        
        
        function createVoigtMatrix(obj)
            theta = obj.Angle;
            u1 = obj.Direction(1);
            u2 = obj.Direction(2);
            u3 = obj.Direction(3);
            
            A(1,1) = ((1 - cos(theta))*u1^2 + cos(theta))^2;
            A(2,2) = ((1 - cos(theta))*u2^2 + cos(theta))^2;
            A(3,3) = ((1 - cos(theta))*u3^2 + cos(theta))^2;
            
            A(4,4) = ((1 - cos(theta))*u2^2 + cos(theta))*((1 - cos(theta))*u3^2 + cos(theta)) - (u1*sin(theta) + u2*u3*(cos(theta) - 1))*(u1*sin(theta) - u2*u3*(cos(theta) - 1));
            A(5,5) = ((1 - cos(theta))*u1^2 + cos(theta))*((1 - cos(theta))*u3^2 + cos(theta)) - (u2*sin(theta) + u1*u3*(cos(theta) - 1))*(u2*sin(theta) - u1*u3*(cos(theta) - 1));
            A(6,6) = ((1 - cos(theta))*u1^2 + cos(theta))*((1 - cos(theta))*u2^2 + cos(theta)) - (u3*sin(theta) + u1*u2*(cos(theta) - 1))*(u3*sin(theta) - u1*u2*(cos(theta) - 1));
            
            A(1,2) = (u3*sin(theta) - u1*u2*(cos(theta) - 1))^2;
            A(2,1) = (u3*sin(theta) + u1*u2*(cos(theta) - 1))^2;
            
            A(1,3) = (u2*sin(theta) + u1*u3*(cos(theta) - 1))^2;
            A(3,1) = (u2*sin(theta) - u1*u3*(cos(theta) - 1))^2;
            
            A(1,4) = -2*(u2*sin(theta) + u1*u3*(cos(theta) - 1))*(u3*sin(theta) - u1*u2*(cos(theta) - 1));
            A(4,1) = -(u3*sin(theta) + u1*u2*(cos(theta) - 1))*(u2*sin(theta) - u1*u3*(cos(theta) - 1));
            
            A(1,5) = -2*(u2*sin(theta) + u1*u3*(cos(theta) - 1))*((1 - cos(theta))*u1^2 + cos(theta));
            A(5,1) = (u2*sin(theta) - u1*u3*(cos(theta) - 1))*((1 - cos(theta))*u1^2 + cos(theta));
            
            A(1,6) = 2*(u3*sin(theta) - u1*u2*(cos(theta) - 1))*((1 - cos(theta))*u1^2 + cos(theta));
            A(6,1) = -(u3*sin(theta) + u1*u2*(cos(theta) - 1))*((1 - cos(theta))*u1^2 + cos(theta));
            
            A(2,3) = (u1*sin(theta) - u2*u3*(cos(theta) - 1))^2;
            A(3,2) = (u1*sin(theta) + u2*u3*(cos(theta) - 1))^2;
            
            A(2,4) = 2*(u1*sin(theta) - u2*u3*(cos(theta) - 1))*((1 - cos(theta))*u2^2 + cos(theta));
            A(4,2) = -(u1*sin(theta) + u2*u3*(cos(theta) - 1))*((1 - cos(theta))*u2^2 + cos(theta));
            
            A(2,5) = -2*(u3*sin(theta) + u1*u2*(cos(theta) - 1))*(u1*sin(theta) - u2*u3*(cos(theta) - 1));
            A(5,2) = -(u1*sin(theta) + u2*u3*(cos(theta) - 1))*(u3*sin(theta) - u1*u2*(cos(theta) - 1));
            
            A(2,6) = -2*(u3*sin(theta) + u1*u2*(cos(theta) - 1))*((1 - cos(theta))*u2^2 + cos(theta));...
            A(6,2) = (u3*sin(theta) - u1*u2*(cos(theta) - 1))*((1 - cos(theta))*u2^2 + cos(theta));
            
            A(3,4) = -2*(u1*sin(theta) + u2*u3*(cos(theta) - 1))*((1 - cos(theta))*u3^2 + cos(theta));
            A(4,3) = (u1*sin(theta) - u2*u3*(cos(theta) - 1))*((1 - cos(theta))*u3^2 + cos(theta));
            
            A(3,5) = 2*(u2*sin(theta) - u1*u3*(cos(theta) - 1))*((1 - cos(theta))*u3^2 + cos(theta));
            A(5,3) = -(u2*sin(theta) + u1*u3*(cos(theta) - 1))*((1 - cos(theta))*u3^2 + cos(theta));
            
            A(3,6) = -2*(u1*sin(theta) + u2*u3*(cos(theta) - 1))*(u2*sin(theta) - u1*u3*(cos(theta) - 1));
            A(6,3) = -(u2*sin(theta) + u1*u3*(cos(theta) - 1))*(u1*sin(theta) - u2*u3*(cos(theta) - 1));
                        
            A(4,5) = (u1*sin(theta) - u2*u3*(cos(theta) - 1))*(u2*sin(theta) - u1*u3*(cos(theta) - 1)) - (u3*sin(theta) + u1*u2*(cos(theta) - 1))*((1 - cos(theta))*u3^2 + cos(theta));
            A(5,4) = (u1*sin(theta) + u2*u3*(cos(theta) - 1))*(u2*sin(theta) + u1*u3*(cos(theta) - 1)) + (u3*sin(theta) - u1*u2*(cos(theta) - 1))*((1 - cos(theta))*u3^2 + cos(theta));
            
            A(4,6) = (u1*sin(theta) + u2*u3*(cos(theta) - 1))*(u3*sin(theta) + u1*u2*(cos(theta) - 1)) + (u2*sin(theta) - u1*u3*(cos(theta) - 1))*((1 - cos(theta))*u2^2 + cos(theta));
            A(6,4) = (u1*sin(theta) - u2*u3*(cos(theta) - 1))*(u3*sin(theta) - u1*u2*(cos(theta) - 1)) - (u2*sin(theta) + u1*u3*(cos(theta) - 1))*((1 - cos(theta))*u2^2 + cos(theta));
            
            A(5,6) = (u2*sin(theta) - u1*u3*(cos(theta) - 1))*(u3*sin(theta) - u1*u2*(cos(theta) - 1)) - (u1*sin(theta) + u2*u3*(cos(theta) - 1))*((1 - cos(theta))*u1^2 + cos(theta));
            A(6,5) = (u2*sin(theta) + u1*u3*(cos(theta) - 1))*(u3*sin(theta) + u1*u2*(cos(theta) - 1)) + (u1*sin(theta) - u2*u3*(cos(theta) - 1))*((1 - cos(theta))*u1^2 + cos(theta));
            
            obj.VoigtMatrix = A;
                
        end
        
        function createVoigtMatrixPlaneStress(obj)
            PSIndex = PlaneStressIndex();
            InPlaneIndex(:,1) = PSIndex.getInPlaneIndex();
            obj.VoigtMatrixPlaneStress = obj.VoigtMatrix(InPlaneIndex,InPlaneIndex);
        end
        
        
    end
      
    
end

