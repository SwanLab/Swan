classdef StressComputer < handle
    
    properties (Access = private)
        C
        dim
        strain
    end
    
    methods (Access = public)
        
        function obj = StressComputer(cParams)
            obj.init(cParams)
        end

        function stress = compute(obj)
            stress = obj.computeStress();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.C      = cParams.C;
            obj.dim    = cParams.dim;
            obj.strain = cParams.strain;
        end

        function stress = computeStress(obj)
            ngaus = obj.dim.ngaus;
            nstre = obj.dim.nstre;
            Cmat  = obj.C;
            strn  = permute(obj.strain,[2 3 1]);
            stress = zeros(size(strn));
            for igaus = 1:ngaus
                for istre=1:nstre
                    for jstre=1:nstre
                        strainIJ = strn(jstre,:,igaus);
                        CmatIJ = Cmat(istre,jstre,:,igaus);
                        Cmatsq = squeeze(CmatIJ);
                        inc = Cmatsq'.*strainIJ;
                        stress(istre,:,igaus) = stress(istre,:,igaus) + inc;
                    end
                end
            end
            stress = permute(stress, [3 1 2]);
        end
        
    end
    
end