classdef DomainFunction < BaseFunction

    properties (GetAccess = public, SetAccess = protected)
        operation
    end
    
    methods (Access = public)
        
        function obj = DomainFunction(cParams)
            obj.init(cParams)
        end
    
    end

   methods (Access = private)

       function init(obj,cParams)
           obj.operation = cParams.operation;
           obj.ndimf = cParams.ndimf;
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

  methods (Access = protected)

        function fxV = evaluateNew(obj, xV)
            fxV = obj.operation(xV);
        end        

    end
   

end