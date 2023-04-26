classdef DesignVariable < handle
    
    properties (Access = public)
        thetavec
    end

    properties (Access = private)
        neuronsPerLayer
        nLayers
    end
    
    methods (Access = public)

        function obj = DesignVariable(s)
            obj.init(s);
            obj.thetavec = obj.computeInitialTheta();
            disp(obj.thetavec);
        end  

    end

    methods (Access = private)

        function init(obj,s)
            obj.neuronsPerLayer = s(2);
            obj.nLayers = length(s(2)); 
        end

        function th = computeInitialTheta(obj)
           nPL    = obj.neuronsPerLayer;
           th     = [];
           for i = 2:obj.nLayers
                if i ~= obj.nLayers
                    b = zeros([1,nPL(i)]) + 0.1;
                else
                    b = zeros([1,nPL(i)]) + 1/nPL(i);
                end
                u = (6/(nPL(i-1)+nPL(i)))^0.5;
                W = (unifrnd(-u,u,[1,nPL(i-1)*nPL(i)]));
                th = [th,W,b];
           end      
       end  
    end
end

