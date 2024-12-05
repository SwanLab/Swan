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
            obj.mesh = cParams.mesh;
        end

   end

   methods (Access = public, Static)

       function f = create(operation, mesh, ndimf)
           if nargin == 3, s.ndimf = ndimf; end
           s.operation = operation;
           s.mesh = mesh;
           f = DomainFunction(s);
       end




   end


    
    


end