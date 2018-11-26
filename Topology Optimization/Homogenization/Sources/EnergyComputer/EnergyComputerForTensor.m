classdef EnergyComputerForTensor < EnergyComputer
    
    properties
    end
    
    
    methods (Access = public)
        
        function obj = EnergyComputerForTensor(C,s)
            obj.init(C,s)
            obj.computeEnergy()
        end
    end
    
    methods (Access = protected)
        
        function computeEnergy(obj)
            C = obj.fourthOrder.getValue();
            s = obj.secondOrder.getValue();
            d = obj.secondOrder.getDimension();
            en = 0;
            for i = 1:d
                for j = 1:d
                    for k = 1:d
                        for l = 1:d
                            sij   = s(i,j);
                            cijkl = C(i,j,k,l);
                            skl   = s(k,l);
                            e = sij*cijkl*skl;
                            en = en + e;
                        end
                    end
                end
            end
            obj.energy = en;
        end
        
    end
end


