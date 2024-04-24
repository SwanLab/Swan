classdef CharacteristicFunctionComputer < handle

    properties (Access = public)
        tfi
        fi
    end

    properties (Access = private)
        levelSet
        t
        p
    end

    methods (Access = public)

        function obj = CharacteristicFunctionComputer(cParams)
            obj.init(cParams)
            obj.compute();
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.p = cParams.p;
            obj.t = cParams.t;
            obj.levelSet = cParams.psi;
        end

        function compute(obj)
            psi = obj.levelSet; 

            n = size(psi,2);
            chi = psi<0;
            
            if n>1
                for i=1:n-1
                    obj.fi(:,i) = (1 - chi(:,i+1)).*prod(chi(:,1:i),2);
                end
                obj.fi(:,n) = prod(chi,2);
            else
                obj.fi(:,1) = chi*1; %multiplied by 1 in order to convert from logical to double
            end
            obj.fi(:,end+1) = (1 - chi(:,1));
            
         %   if nargout==2
                for i=1:n
                [tXi(i,:),~] = integ_exact(obj.t,obj.p,psi(:,i));
                end
                tXi = 1 - tXi;
                if n>1
                    for i=1:n-1
                        obj.tfi(i,:) = (1 - tXi(i+1,:)).*prod(tXi(1:i,:),1);
                    end
                    obj.tfi(n,:) = prod(tXi,1);
                else
                    obj.tfi(1,:) = tXi;
                end
                obj.tfi(end+1,:) = (1 - tXi(1,:));    
         %   end
                    end
                end

    
end