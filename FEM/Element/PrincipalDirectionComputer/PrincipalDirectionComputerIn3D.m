classdef PrincipalDirectionComputerIn3D < PrincipalDirectionComputer
    
    methods (Access = public, Static)
        
        function obj = PrincipalDirectionComputerIn3D(cParams)
            obj.ndim = 3;
            obj.init()
        end
        
    end
    
    methods (Access = public)
        
        function compute(obj,tensor)
            s1  = squeeze(tensor(1,1,15));
            s2  = squeeze(tensor(1,2,15));
            s3  = squeeze(tensor(1,3,15));
            s12 = squeeze(tensor(1,4,15));
            s13 = squeeze(tensor(1,5,15));
            s23 = squeeze(tensor(1,6,15));
            for i = 1:obj.ndim
                for j = 1:obj.ndim
                    d = obj.directionFunction{i,j}(s1,s2,s3,s12,s13,s23);
                    obj.direction(i,j,:) = real(d);                    
                    obj.assertSmallImaginaryValue(d);
                end
            end
        end
        
    end
    
    methods (Access = protected)
        
        function obtainEigenVectors(obj)
            s1 = sym('s1','real');
            s2 = sym('s2','real');
            s3 = sym('s3','real');
            s12 = sym('s12','real');
            s13 = sym('s13','real');
            s23 = sym('s23','real');
            S = [s1 s12 s13; s12 s2 s23; s13 s23 s3];
            [vS,~] = eig(S);
            obj.eigenVectors = vS;
        end
                
    end
    
    methods (Access = private, Static)
        
        function assertSmallImaginaryValue(d)
            imagD = imag(d);
            assert(norm(imagD(:)) < 1e-1)
        end        
        
    end
    
end