classdef InitialLevelSetComputer < handle

    properties (Access = public)
        psi
    end

    properties (Access = private)
        E
        p
    end

    methods (Access = public)

        function obj = InitialLevelSetComputer(mat,cParams)
            obj.init(mat,cParams);
            obj.computeInitialGuess();
        end
        
    end

    methods (Access = private)

        function init(obj,mat,cParams)
            obj.E(1) = mat.A.young;
            obj.E(2) = mat.B.young;
            obj.E(3) = mat.C.young;
            obj.E(4) = mat.D.young;

            obj.p = cParams.p; 
        end

        function computeInitialGuess(obj)
            obj.psi = ones(size(obj.p,2),length(obj.E)-1); 
            obj.psi(:,1) = -1;
        end
    end
end