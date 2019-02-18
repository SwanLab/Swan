classdef FeasibleIndexComputer < handle
    
    properties (SetAccess = private, GetAccess = public) 
        index
    end
    
    properties (Access = private)
        mx
        my
        nmx
        nmy
        chi
        rho
        iIndex
        jIndex         
    end
    
    
    methods (Access = public)
        
        function obj = FeasibleIndexComputer(d)
            obj.init(d);
            obj.computeFeasibleIJindex();
            obj.computeGlobalIndex();
        end
    end
    
    
    methods (Access = private)
        
        function init(obj,d)
            obj.mx = d.mx;            
            obj.my = d.my;            
            obj.chi = d.chi;
            obj.rho = d.rho;
            obj.nmx = length(obj.mx);
            obj.nmy = length(obj.my);
        end
        
          function computeFeasibleIJindex(obj)
            k = 0;
            for i = 1:obj.nmx
                for j = 1:obj.nmy
                    if obj.isFeasible(i,j)
                        k = k+1;
                        obj.iIndex(k) = i;
                        obj.jIndex(k) = j;
                    end
                end
            end
        end
        
        function itIs = isFeasible(obj,i,j)
            itIs = obj.isFirstCondSatisfied(i,j) &&...
                obj.isSecondCondSatisfied(i,j) && ...
                obj.isThirdCondSatisfied(i,j);
        end
        
        function itIs = isFirstCondSatisfied(obj,i,j)
            itIs = obj.chi(i,j) <= 1;
        end
        
        function itIs = isSecondCondSatisfied(obj,i,j)
            itIs = 1 - obj.rho(i,j) <= obj.chi(i,j);
        end
        
        function itIs = isThirdCondSatisfied(obj,i,j)
            itIs = obj.chi(i,j) <= 1/(1-obj.rho(i,j));
        end 
        
        function computeGlobalIndex(obj)
            i = obj.iIndex();
            j = obj.jIndex();
            obj.index = sub2ind([obj.nmx,obj.nmy],i,j);            
        end
        
    end
    
end