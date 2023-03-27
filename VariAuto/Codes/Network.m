classdef Network < handle
 
    properties (GetAccess = public, SetAccess = private)
       data
       thetavec
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

       function self = Network(varargin)
           self.init(varargin);
           self.computeInitialTheta();
       end

       function computeInitialTheta(self)
           nPL    = self.neuronsPerLayer;
           th     = [];
           for i = 2:self.nLayers
                if i ~= self.nLayers
                    b = zeros([1,nPL(i)]) + 0.1;
                else
                    b = zeros([1,nPL(i)]) + 1/nPL(i);
                end
                u = (6/(nPL(i-1)+nPL(i)))^0.5;
                W = (unifrnd(-u,u,[1,nPL(i-1)*nPL(i)]));
                th = [th,W,b];
           end      
           self.thetavec = th;
       end
   end

   methods (Access = private)

       function init(self,s)
           self.data = s{1};
           self.neuronsPerLayer = s{2};
           self.nLayers = length(s{2});
           if length(s) <= 2
               self.Costtype = '-loglikelihood';
               self.HUtype = 'ReLU';
               self.OUtype = 'softmax';
               self.lambda = 0;
           else
               self.Costtype = s{3};
               self.HUtype = s{4};
               self.OUtype = s{5};
               self.lambda = s{6};
           end
       end  
   end   

   methods 
       function value = get.layer(self)
            nPL = self.neuronsPerLayer;
            last = 1;
            value = cell(self.nLayers-1,1);
            for i = 2:self.nLayers
                aux = nPL(i)*nPL(i-1) + nPL(i);
                next = last + aux;
                theta_i = self.thetavec(last:next-1);
                value{i-1} = Layer(theta_i,nPL(i-1),nPL(i));
                last = next;
            end
       end
   end
end