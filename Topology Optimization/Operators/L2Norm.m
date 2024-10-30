classdef L2Norm < handle


    methods (Access = public, Static)
        
        function sp = compute(m,f)
            sp = Integrator.compute(f.*f,m,2);
            sp = sqrt(sp);
        end

    end

end