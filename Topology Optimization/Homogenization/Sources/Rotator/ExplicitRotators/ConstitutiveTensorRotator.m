classdef ConstitutiveTensorRotator < handle
    
    properties (Access = private)
        prodFunc
        Rot
        Cs
    end
    
    
    methods (Access = public)
        
        function obj = ConstitutiveTensorRotator()
            obj.createSymbolicRotator();
            obj.createSymbolicConstitutiveTensor();
            obj.createProductFunction();
        end
        
        function c = rotate(obj,C,R)
            R11 = squeeze(R(1,1,:));
            R12 = squeeze(R(1,2,:));
            R13 = squeeze(R(1,3,:));
            R21 = squeeze(R(2,1,:));            
            R22 = squeeze(R(2,2,:));
            R23 = squeeze(R(2,3,:));            
            R31 = squeeze(R(3,1,:));
            R32 = squeeze(R(3,2,:));
            R33 = squeeze(R(3,3,:));
            
            C11 = squeeze(C(1,1,:));
            C12 = squeeze(C(1,2,:));
            C13 = squeeze(C(1,3,:));
            C22 = squeeze(C(2,2,:));
            C23 = squeeze(C(2,3,:));
            C33 = squeeze(C(3,3,:));
            
            c(1,1,:) = obj.prodFunc{1,1}(C11,C12,C13,C22,C23,C33,R11,R21,R31);
            c(1,2,:) = obj.prodFunc{1,2}(C11,C12,C13,C22,C23,C33,R11,R12,R21,R22,R31,R32);
            c(1,3,:) = obj.prodFunc{1,3}(C11,C12,C13,C22,C23,C33,R11,R13,R21,R23,R31,R33);
            c(2,2,:) = obj.prodFunc{2,2}(C11,C12,C13,C22,C23,C33,R12,R22,R32);
            c(2,3,:) = obj.prodFunc{2,3}(C11,C12,C13,C22,C23,C33,R12,R13,R22,R23,R32,R33);
            c(3,3,:) = obj.prodFunc{3,3}(C11,C12,C13,C22,C23,C33,R13,R23,R33);
            c(2,1,:) = c(1,2,:);
            c(3,1,:) = c(1,3,:);    
            c(3,2,:) = c(2,3,:); 
        end
        
    end
    
    methods (Access = private)
        
        function createSymbolicRotator(obj)
            R11 = sym('R11','real');
            R12 = sym('R12','real');
            R13 = sym('R13','real');
            R21 = sym('R21','real');
            R22 = sym('R22','real');
            R23 = sym('R23','real');
            R31 = sym('R31','real');
            R32 = sym('R32','real');
            R33 = sym('R33','real');            
            obj.Rot = [R11 R12 R13;
                      R21 R22 R23;
                      R31 R32 R33];            
        end
        
        function createSymbolicConstitutiveTensor(obj)
            C11 = sym('C11','real');
            C12 = sym('C12','real');
            C13 = sym('C13','real');
            C22 = sym('C22','real');
            C23 = sym('C23','real');
            C33 = sym('C33','real');            
            obj.Cs = [C11 C12 C13;
                C12 C22 C23;
                C13 C23 C33];            
        end
        

        function createProductFunction(obj)            
            R = obj.Rot;
            C = obj.Cs;
            prod = simplify(R'*C*R);
            for i = 1:3
                for j=1:3                    
                   obj.prodFunc{i,j} = matlabFunction(prod(i,j));
                end
            end
        end
        
    end
    
end