classdef PhaseFieldDegradationInterpolator < handle

    methods (Static, Access = public)

        function [eq,deq,d2eq] = compute(type,f0)
            syms x;
            switch type
                case 'AT'
                    deg = ((1-x).^2);
            end
            
            eq   = matlabFunction(deg,'Vars',x);
            eq   = @(phi) eq(phi).*f0;
            deq  = matlabFunction(diff(deg),'Vars',x);
            deq  = @(phi) deq(phi).*f0;
            d2eq = matlabFunction(diff(diff(deg)),'Vars',x);
            d2eq = @(phi) d2eq(phi).*f0;
        end

    end
end