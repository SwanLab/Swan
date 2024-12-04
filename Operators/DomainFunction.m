classdef DomainFunction < BaseFunction

    properties (GetAccess = public, SetAccess = protected)
        operation
    end
    
    methods (Access = public)
        
        function obj = DomainFunction(cParams)
            obj.init(cParams)
        end
        
        function r = evaluate(obj,xV)
            r = obj.operation(xV);
        end

    end

   methods (Access = private)

        function init(obj,cParams)
            obj.operation = cParams.operation;
            if isfield(cParams,'ndimf')
                obj.ndimf = cParams.ndimf;
            else
                obj.ndimf = 1;
            end

            if isfield(cParams,'mesh')
                obj.mesh = cParams.mesh;
            else
                obj.mesh = [];
            end
        end

    end
    
    


end