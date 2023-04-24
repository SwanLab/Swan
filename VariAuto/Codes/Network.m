classdef Network < handle
 
    properties (Access = public)
        thetavec
    end

    properties (GetAccess = public, SetAccess = private)
       data       
       neuronsPerLayer
       nLayers
       lambda
       Costtype
       HUtype
       OUtype
    end

    properties (Dependent)
        layer
    end

   methods (Access = public)

       function obj = Network(varargin)
           obj.init(varargin);
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

   methods (Access = private)

       function init(obj,s)
           obj.data = s{1};
           obj.neuronsPerLayer = s{2};
           obj.nLayers = length(s{2});
           if length(s) <= 2
               obj.Costtype = '-loglikelihood';
               obj.HUtype = 'ReLU';
               obj.OUtype = 'softmax';
               obj.lambda = 0;
           else
               obj.Costtype = s{3};
               obj.HUtype = s{4};
               obj.OUtype = s{5};
               obj.lambda = s{6};
           end
       end  
   end   

   methods 
       function value = getLayer(obj)
            nPL = obj.neuronsPerLayer;
            last = 1;
            value = cell(obj.nLayers-1,1);
            for i = 2:obj.nLayers
                aux = nPL(i)*nPL(i-1) + nPL(i);
                next = last + aux;
                theta_i = obj.thetavec(last:next-1);
                value{i-1} = Layer(theta_i,nPL(i-1),nPL(i));
                last = next;
            end
       end
   end
end