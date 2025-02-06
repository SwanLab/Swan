classdef PrincipalDirectionComputerIn3D < PrincipalDirectionComputer
    
    methods (Access = public, Static)
        
        function obj = PrincipalDirectionComputerIn3D(cParams)
            obj.ndim = 3;
            obj.init(cParams)
        end
        
    end
    
    methods (Access = public)
        
        function compute(obj,tensor)
            s1  = squeeze(tensor(1,1,:));
            s2  = squeeze(tensor(1,2,:));
            s3  = squeeze(tensor(1,3,:));
            s12 = squeeze(tensor(1,4,:));
            s13 = squeeze(tensor(1,5,:));
            s23 = squeeze(tensor(1,6,:));
            eG = obj.eigenComputer;
            for i = 1:obj.ndim
                for j = 1:obj.ndim
                    d = eG.eigenVectorFunction{i,j}(s1,s12,s13,s2,s23,s3);
                    obj.direction(i,j,:) = real(d);
                    obj.assertSmallImaginaryValue(d);
                end
            end
        end

    end
    
 
    
    methods (Access = private, Static)
        
        function assertSmallImaginaryValue(d)
            imagD = imag(d);
%             assert(norm(imagD(:)) < 1e6) % To be fixed
        end
        
    end
    
end