classdef ElasticProblemDisp < ElasticProblem

    properties (Access = public)
        sizePer
    end

    methods (Access = public)

        function Ch = computeStressHomog(obj)
            nstre = obj.material.nstre;
            basis = diag(ones(nstre,1));
            Ch    = zeros(nstre,nstre);

            for istre = 1:nstre
                obj.vstrain = basis(istre,:);
                obj.boundaryConditions.updateBC(obj.vstrain);
                obj.solve();
                perDOFslave = obj.boundaryConditions.periodic_constrained;
                obj.sizePer = size(perDOFslave, 1);
                nEqperType  = obj.sizePer/4;
                L           = obj.variables.LangMult;
                Lx          = 0;
                Ly          = 0;
                Lxy         = 0;
                d1          = nEqperType*3;
                LPer        = L(1:d1);
                LVer        = L(d1+1:end);
                for i = 1:nEqperType
                    Lx  = Lx + LPer(i);
                end
                for i = nEqperType+1:2*nEqperType
                    Lxy = Lxy + LPer(i);
                end
                for i = 2*nEqperType+1:3*nEqperType
                    Ly  = Ly + LPer(i);
                end
    
                if obj.vstrain(1) == 1
                    Lx  = Lx + LVer(1) + LVer(2) + LVer(3) + LVer(5);
                    Ly  = Ly + LVer(4) + LVer(7);
                elseif obj.vstrain(2) == 1
                    Ly  = Ly + LVer(1) + LVer(2) + LVer(4) + LVer(6);
                    Lx  = Lx + LVer(3) + LVer(7);                
                else
                    Lxy = Lxy + LVer(1) + LVer(2);
                end
                Ch(istre, :) = [Lx; Ly; Lxy];                
            end
            obj.variables.Chomog  = -Ch;

        end
    end

end