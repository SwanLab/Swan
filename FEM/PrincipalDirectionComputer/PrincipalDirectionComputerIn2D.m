classdef PrincipalDirectionComputerIn2D < PrincipalDirectionComputer
    
    properties (Access = private)
       avarageTensor               
    end
    
    methods (Access = public)
        
        function obj = PrincipalDirectionComputerIn2D(cParams)
            obj.ndim = 2;
            obj.init(cParams);
        end
        
    end
    
    methods (Access = public)
        
        function compute(obj,tensor)
            obj.computeAvarageTensor(tensor);
            s1  = squeeze(obj.avarageTensor(1,1,:));
            s2  = squeeze(obj.avarageTensor(1,2,:));
            s12 = squeeze(obj.avarageTensor(1,3,:));
            eG = obj.eigenComputer;
            for i = 1:2
                for j = 1:2
                    obj.direction(i,j,:) = eG.eigenVectorFunction{i,j}(s1,s12,s2);
                end
                obj.principalStress(i,:) = eG.eigenValueFunction{i}(s1,s12,s2);
            end
        end
        
    end
    
    
    methods (Access = private)
        
        function computeAvarageTensor(obj,tensor)
            t = zeros(1,size(tensor,2),size(tensor,3));
            ngaus = size(tensor,1);
            for igaus = 1 : ngaus
               t(1,:,:) = t(1,:,:) + tensor(igaus,:,:)/ngaus;
            end
            obj.avarageTensor = t;
        end
        
    end
    
end