classdef DesVar < handle
    
    properties (Access = public)
        thetavec
    end

    properties (Access = protected)
        network
    end
    
    methods (Access = public)
        function obj = DesVar(varargin)
            obj.init(varargin);
        end

        function computeInitialTheta(obj)
           nPL    = obj.network.neuronsPerLayer;
           th     = [];
           for i = 2:obj.network.nLayers
                if i ~= obj.network.nLayers
                    b = zeros([1,nPL(i)]) + 0.1;
                else
                    b = zeros([1,nPL(i)]) + 1/nPL(i);
                end
                u = (6/(nPL(i-1)+nPL(i)))^0.5;
                W = (unifrnd(-u,u,[1,nPL(i-1)*nPL(i)]));
                th = [th,W,b];
           end      
           obj.thetavec = th;
       end
    end

    methods (Access = private)
        function init(obj,s)
           obj.network = s{1};
       end  
    end
end

