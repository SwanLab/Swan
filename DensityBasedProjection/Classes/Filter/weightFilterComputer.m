classdef weightFilterComputer < handle
    properties (Access = public)
        H
        Hs        
    end
    properties (Access = private)
        Xnumber
        Ynumber
        minimunInfluenceRadios
    end

    methods (Access = public)
        function obj = weightFilterComputer(cParams)
            obj.inputData(cParams);
        end
        function compute(obj)
            obj.computeWeight();
        end
    end
    methods (Access = private)
        function obj = inputData(obj,cParams)
            obj.Xnumber =  cParams.Xnumber;
            obj.Ynumber =  cParams.Ynumber;
            obj.minimunInfluenceRadios = cParams.minimunInfluenceRadios;
        end
        function computeWeight(obj)
            iH = ones(obj.Xnumber*obj.Ynumber*(2*(ceil(obj.minimunInfluenceRadios)-1)+1)^2,1);
            jH = ones(size(iH));
            sH = zeros(size(iH));
            k  = 0;
            for i1 = 1:obj.Xnumber
                for j1 = 1:obj.Ynumber
                    e1 = (i1-1)*obj.Ynumber+j1;
                    for i2 = max(i1-(ceil(obj.minimunInfluenceRadios)-1),1):min(i1+(ceil(obj.minimunInfluenceRadios)-1),obj.Xnumber)
                        for j2 = max(j1-(ceil(obj.minimunInfluenceRadios)-1),1):min(j1+(ceil(obj.minimunInfluenceRadios)-1),obj.Ynumber)
                            e2 = (i2-1)*obj.Ynumber+j2;
                            k = k+1;
                            iH(k) = e1;
                            jH(k) = e2;
                            sH(k) = max(0,obj.minimunInfluenceRadios-sqrt((i1-i2)^2+(j1-j2)^2));
                        end
                    end
                end
            end
            obj.H        = sparse(iH,jH,sH);
            obj.Hs       = sum(obj.H,2);            
        end 
    end
end
