classdef SectionVariablesComputer < handle
    
    properties (GetAccess = public, SetAccess = protected)
        designVariable
        mesh
        nDesVarElem
    end
    
    properties (Access = private)
 
    end
    
    methods (Access = public, Static)
        
        function section = create(cParams)
            f = SectionVariableFactory();
            section = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
        
        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
            obj.mesh = cParams.mesh;
            obj.nDesVarElem = cParams.designVariable.nDesignVar;             
        end

        function value = getSingleValue(obj)
            x = obj.designVariable.value;
            N = obj.mesh.nelem;
            value = x(1:N);
        end

        function [a,b] = getDoubleValue(obj)
            x = obj.designVariable.value;
            N = obj.mesh.nelem;
            a = x(1:N);
            b = x(N+1:2*N);
        end
        
        function [a,b,c,d] = getFourValues(obj)
            x = obj.designVariable.value;
            N = obj.mesh.nelem;
            a = x(1:N);
            b = x(N+1:2*N);
            c = x(2*N+1:3*N);
            d = x(3*N+1:4*N);
        end
    end
    
end