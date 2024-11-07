classdef Norm < handle

  methods (Access = public, Static)
        
        function sp = computeL2(m,f)
           sp = ScalarProduct.computeL2(m,f,f);
           sp = sqrt(sp);
        end

        function sp = computeH1(m,f,eps)
            sp = ScalarProduct.computeH1(m,f,f,eps);
        end

    end

end