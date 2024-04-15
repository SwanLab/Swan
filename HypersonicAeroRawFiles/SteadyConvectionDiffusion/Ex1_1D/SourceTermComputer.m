classdef SourceTermComputer < handle

    properties (Access = private)
        fHandle
    end

    methods (Static, Access = public)
        function obj = create(fHandle)
            s.source = fHandle;
            obj      = SourceTermComputer(s);
        end
    end

    methods (Access = public)
        function obj = SourceTermComputer(cParams)
            obj.init(cParams);
        end

        function s = compute(obj,x)
            s = obj.fHandle(x);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.fHandle = cParams.source;
        end
    end 
end