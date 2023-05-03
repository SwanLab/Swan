classdef Network < handle
 
    properties (Access = public)
        theta
    end

    properties (GetAccess = public, SetAccess = private)
       neuronsPerLayer
       nLayers
       lambda
       Costtype
       HUtype
       OUtype
    end

    properties (Access = private)
       data  
    end

    properties (Dependent)
        layer
        W
        b
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
                    getB = zeros([1,nPL(i)]) + 0.1;
                else
                    getB = zeros([1,nPL(i)]) + 1/nPL(i);
                end
                u = (6/(nPL(i-1)+nPL(i)))^0.5;
                getW = (unifrnd(-u,u,[1,nPL(i-1)*nPL(i)]));
                th = [th,getW,getB];
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

       function s = getPrev(obj,thv,prevL,nextL)
            s.theta = thv;
            s.prevL = prevL;
            s.nextL = nextL;
            s.b     = obj.getB(thv,prevL,nextL);
            s.W     = obj.getW(thv,prevL,nextL);
        end
   end   

   methods 
       function value = getLayer(obj,thetavec)
            nPL = obj.neuronsPerLayer;
            last = 1;
            value = cell(obj.nLayers-1,1);
            for i = 2:obj.nLayers
                aux = nPL(i)*nPL(i-1) + nPL(i);
                next = last + aux;
                theta_i = thetavec(last:next-1);                
                value{i-1} = obj.getPrev(theta_i,nPL(i-1),nPL(i));
                last = next;
            end
       end

       function value = getW(obj,thv,prevL,nextL)
            aux = thv(1:prevL*nextL);
            value = reshape(aux,[prevL,nextL]);
       end

        function value = getB(obj,thv,prevL,nextL)
            value = thv(prevL*nextL+1:end);
        end
   end

end