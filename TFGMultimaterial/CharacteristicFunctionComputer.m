classdef CharacteristicFunctionComputer < handle

    properties (Access = private)
        t
        p
        designVariable
        charfuncTfi
        charfuncFi
        mesh
    end

    methods (Access = public)

        function obj = CharacteristicFunctionComputer(cParams)
            obj.init(cParams);
            obj.computeDomainFunctionTfi();
            obj.computeDomainFunctionFi();
            obj.computeFiandTfi();
        end


        function [fi, tfi] = computeFiandTfi(obj)
            cParams.chi = obj.charfuncTfi;
            tfi = computeFunction(cParams,obj);
            
            cParams.chi = obj.charfuncFi;
            fi = computeFunction(cParams,obj);     
        end
    end
    

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.m;
            obj.designVariable = cParams.designVariable;
        end

        function x = computeFunction(cParams,obj)
            n = 3; % size psi
            chi = cParams.chi;
            if n>1
                for i=1:n-1
                    phi(:,i) = (1 - chi(:,i+1)).*prod(chi(:,1:i),2);
                end
                phi(:,n) = prod(chi,2);
            else
                phi(:,1) = chi*1; %multiplied by 1 in order to convert from logical to double
            end
            phi(:,end+1) = (1 - chi(:,1));

            s.type = 'Full';
            s.mesh = obj.mesh;
            s.plotting = false;
            s.fValues = phi;
            s.order   = 'P1';

            fValues   = LagrangianFunction(s); 
            x = fValues.fValues;
            x = x';
        end

        function computeDomainFunctionTfi(obj)
            x1 = obj.designVariable.designVariable{1};
            chi1 = x1.obtainDomainFunction();
            chi1 = chi1.project('P0');

            x2 = obj.designVariable.designVariable{2};
            chi2 = x2.obtainDomainFunction();
            chi2 = chi2.project('P0');

            x3 = obj.designVariable.designVariable{3};
            chi3 = x3.obtainDomainFunction();
            chi3 = chi3.project('P0');

            obj.charfuncTfi(:,1) = chi1.fValues;
            obj.charfuncTfi(:,2) = chi2.fValues;
            obj.charfuncTfi(:,3) = chi3.fValues;
        end

        function computeDomainFunctionFi(obj)
            x1 = obj.designVariable.designVariable{1};
            chi1 = x1.obtainDomainFunction();
            chi1 = chi1.project('P1');

            x2 = obj.designVariable.designVariable{2};
            chi2 = x2.obtainDomainFunction();
            chi2 = chi2.project('P1');

            x3 = obj.designVariable.designVariable{3};
            chi3 = x3.obtainDomainFunction();
            chi3 = chi3.project('P1');

            obj.charfuncFi(:,1) = chi1.fValues;
            obj.charfuncFi(:,2) = chi2.fValues;
            obj.charfuncFi(:,3) = chi3.fValues;
        end

    end

end