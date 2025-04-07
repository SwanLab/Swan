classdef PrincipalDirectionComputerIn3D < PrincipalDirectionComputer
    
    methods (Access = public, Static)
        
        function obj = PrincipalDirectionComputerIn3D(cParams)
            obj.ndim = 3;
            obj.init(cParams)
        end
        
    end
    
    methods (Access = public)
        
        function [dF,pF] = compute(obj,tensor)
            s = tensor.fValues;            
            s1  = squeeze(s(:,1));
            s2  = squeeze(s(:,2));
            s3  = squeeze(s(:,3));
            s12 = squeeze(s(:,4));
            s13 = squeeze(s(:,5));
            s23 = squeeze(s(:,6));
            eG = obj.eigenComputer;
            for i = 1:obj.ndim
                for j = 1:obj.ndim
                    d(i,j,:)  = eG.eigenVectorFunction{i,j}(s1,s12,s13,s2,s23,s3);
                end
                p(i,:) = eG.eigenValueFunction{i}(s1,s12,s13,s2,s23,s3);                
            end
            pF = obj.createP1Function(p',tensor.mesh);
            for j = 1:2
                dF{j} = obj.createP1Function(squeeze(d(:,j,:))',tensor.mesh);
            end
            
        end

    end
    
 
    
    methods (Access = private)

        function f = createP1Function(obj,fV,mesh)
            s.fValues = fV;
            s.mesh    = mesh;
            s.order   = 'P1';
            f = LagrangianFunction(s);
        end        
        
    end
    
end