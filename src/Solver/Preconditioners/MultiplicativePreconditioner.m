classdef MultiplicativePreconditioner < handle
    
    properties (Access = private)
        PreconditionerList
    end
    
    methods (Access = public)
        
        function obj = MultiplicativePreconditioner(cParams)
            obj.init(cParams)
        end
        
        function z = apply(obj,A,r)
            z = 0;
            for k = 1:length(obj.PreconditionerList)
                Pk = obj.PreconditionerList{k};
                zk = Pk.apply(A,r);
                r  = r - A*(zk);
                z  = z + zk;
            end
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.PreconditionerList = cParams.varargin;            
        end
        
    end
    
end