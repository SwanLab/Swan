classdef FeFunction < handle
    % NOTE
    % Go to P1Function.m for The Class Formerly Known As FeFunction
    % Eventually should extend Function/Field or something like that, to
    % account for other types of functions (eg. L2)
    
    properties (Constant, Access = public)
        fType = 'FE'
    end

    properties (Access = public)
       ndimf
       order
       fValues
    end

    properties (Access = private)
        
    end
    
    methods (Access = public)

        function obj = FeFunction()
        end

    end

    methods (Static, Access = public)

        function obj = create(cParams)
            fun = FunctionFactory();
            obj = fun.create(cParams);
        end
        
    end

    methods (Access = private)
        
        function init(obj, cParams)

        end

    end

end

