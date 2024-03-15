classdef DesignVariable < handle
    
    properties (GetAccess = public, SetAccess = protected)
        fun
        type        
    end

    properties (Access = protected)
        mesh
        nVariables
        isFixed
    end
    
    properties (Access = private)
        funOld
    end
    
    methods (Access = public, Static)
        
        function designVariable = create(cParams)
            f = DesignVariableFactory();
            designVariable = f.create(cParams);
        end

    end

    methods (Access = public)

        function x = getValue(obj)
            nVar = obj.nVariables;
            nx = length(obj.fun{1});
            x = zeros(nVar*nx,1);
            for ivar = 1:nVar
                i0 = nx*(ivar-1) + 1;
                iF = nx*ivar;
                xs = obj.fun{ivar};
                x(i0:iF) = xs;
            end
        end

        function ndofs = getDofs(obj)
            nVar = obj.nVariables;    
            ndofs = 0;
            for ivar = 1:nVar
                xs = obj.fun{ivar};
                ndofs = ndofs + xs.nDofs;
            end
        end
        
        function updateOld(obj)
            nVar = obj.nVariables;                
            for ivar = 1:nVar
                obj.funOld{ivar} = obj.fun{ivar}.copy();
            end                        
        end
        
        function norm = computeL2normIncrement(obj)
           incFun = obj.fun-obj.funOld;
           nIncX  = Norm.computeL2(obj.mesh,incFun);
           nX0    = Norm.computeL2(obj.mesh,obj.funOld);
           norm   = nIncX/nX0;
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.type = cParams.type;
            obj.mesh = cParams.mesh;
            obj.fun  = cParams.fun;
            if isfield(cParams,'isFixed')
              obj.isFixed = cParams.isFixed;
            end
        end

    end

end