classdef MicroDamageParams < DesignVariable
    
    methods (Access = public)
        
        function obj = MicroDamageParams(cParams)
            obj.nVariables = 1;
            obj.init(cParams)
        end
        
        function update(obj,x)
            x  = obj.splitDesignVariable(x);
            for ivar = 1:obj.nVariables                
               obj.fun{ivar}         = obj.fun{ivar}.copy();
               obj.fun{ivar}.fValues = x{ivar};
            end
        end
        
    end
    
    methods (Access = private)
        
        function xS = splitDesignVariable(obj,x)
            nVar = obj.nVariables;
            nx = length(x)/nVar;
            xS = cell(nVar,1);
            for ivar = 1:nVar
                i0 = nx*(ivar-1) + 1;
                iF = nx*ivar;
                xs = x(i0:iF);
                xS{ivar} = xs;
            end
        end

    end
    
end