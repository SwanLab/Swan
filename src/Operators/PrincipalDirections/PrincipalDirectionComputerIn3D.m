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
            s1  = squeeze(s(1:6:end));
            s2  = squeeze(s(2:6:end));
            s3  = squeeze(s(3:6:end));
            s12 = squeeze(s(4:6:end));
            s13 = squeeze(s(5:6:end));
            s23 = squeeze(s(6:6:end));
            eG = obj.eigenComputer;
            for i = 1:obj.ndim
                for j = 1:obj.ndim
                    d(i,j,:)  = eG.eigenVectorFunction{i,j}(s1,s12,s13,s2,s23,s3);
                end
                p(i,:) = eG.eigenValueFunction{i}(s1,s12,s13,s2,s23,s3);                
            end
            pV = reshape(p,[],1);
            pF = obj.createP1Function(pV,tensor);
            for j = 1:2
                dV = reshape(squeeze(d(:,j,:)),[],1);
                dF{j} = obj.createP1Function(dV,tensor);
            end
            
        end

    end
    
    methods (Access = private)

        function f = createP1Function(~,fV,fun)
            f = LagrangianFunction.create(fun.mesh,fun.ndimf,'P1');
            f.setFValues(fV);
        end        
        
    end
    
end