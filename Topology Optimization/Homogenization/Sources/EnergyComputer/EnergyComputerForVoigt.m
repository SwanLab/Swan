classdef EnergyComputerForVoigt < EnergyComputer
    
    properties
    end
    
    methods (Access = public)
        
        function obj = EnergyComputerForVoigt(C,s)
            obj.init(C,s)
            obj.computeEnergy()
        end
    end
    
    methods (Access = protected)
        
        function computeEnergy(obj)
            C = obj.fourthOrder.getValue();
            s = obj.secondOrder.getValue();
            d = obj.secondOrder.getVoigtDimension();
            en = 0;
            for i = 1:d
                for j = 1:d
                    si  = s(i);
                    Cij = C(i,j);
                    sj = s(j);
                    e = si*Cij*sj;
                    en = en + e;
                end
            end
            obj.energy = en;
        end
        
    end
end