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

        function [dF,pF] = compute(obj,tensor)
            s = tensor.fValues;
            s1  = squeeze(s(1:3:end));
            s2  = squeeze(s(2:3:end));
            s12 = squeeze(s(3:3:end));
            eG = obj.eigenComputer;
            for i = 1:2
                for j = 1:2
                    d(i,j,:) = eG.eigenVectorFunction{i,j}(s1,s12,s2);
                end
                p(i,:) = eG.eigenValueFunction{i}(s1,s12,s2);
            end
            pV = reshape(p,[],1);
            pF = obj.createP1Function(pV,tensor);
            for j = 1:2
                dV = reshape(squeeze(d(:,j,:)),[],1);
                dF{j} = obj.createP1Function(dV,tensor);
            end
        end
        
        function [d,p] = computeFromP0(obj,tensor)
            tensor = obj.computeAverageTensor(tensor);
            [d,p] = obj.compute(tensor);
        end

        
    end
    
    
    methods (Access = private)

        function f = createP1Function(~,fV,fun)
            f = LagrangianFunction.create(fun.mesh,fun.ndimf,'P1');
            f.setFValues(fV);
        end

        function s = computeAverageTensor(obj,tensor)            
            t = zeros(1,size(tensor,2),size(tensor,3));
            ngaus = size(tensor,1);
            for igaus = 1 : ngaus
               t(1,:,:) = t(1,:,:) + tensor(igaus,:,:)/ngaus;
            end
            s(:,1) = squeeze(t(1,1,:));
            s(:,2) = squeeze(t(1,2,:));
            s(:,3) = squeeze(t(1,3,:));            
        end
        
    end
    
end