classdef L2Norm < handle


    methods (Access = public, Static)
        
        function sp = compute(m,f)
            q.mesh         = m;
            q.quadType     = 'QUADRATIC';
            q.type         = 'ScalarProduct';
            int            = Integrator.create(q);
            sp = int.compute(f,f);
        end

    end

end