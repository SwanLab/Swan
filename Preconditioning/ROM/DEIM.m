classdef DEIM < handle
    
    properties (GetAccess = public, SetAccess = private)
        basis
        indices
    end
    
    properties (Access = private)
        data
    end
    
    properties (Access = private)
        threshold
        nBasis
    end
    
    methods (Access = public)
        
        function obj = DEIM(data)
            obj.init(data)
            obj.computeBasis();
            obj.computeMagicPoints();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,data)
            obj.data = data;
            obj.threshold = 0.999;
        end
        
        function computeBasis(obj)
            [U,S,~]    = svd(obj.data,"econ");
            Svec       = diag(S);
            total      = sum(Svec.^2);
            recovered  = 0;
            i          = 1;
            while recovered < obj.threshold
                recovered = recovered + Svec(i)^2/total;
                n=i;
                i=i+1;
            end
            obj.nBasis = n;
            obj.basis = U(:,1:obj.nBasis);
        end
        
        function computeMagicPoints(obj)
            [~, idx(1)] = max(obj.basis(:,1));
            for i = 2:obj.nBasis
                U_sub = obj.basis(idx(1:i-1), 1:i-1);
                c = U_sub \ obj.basis(idx(1:i-1), i);
                r = obj.basis(:,i) - obj.basis(:,1:i-1) * c;
                [~, idx(i)] = max(abs(r));               
            end
            obj.indices = idx;
        end
        
    end
    
end