classdef MultimaterialGradientComputer < handle

    properties (Access = private)
        mesh
        designVariable
    end

    methods (Access = public)
        function obj = MultimaterialGradientComputer(cParams)
            obj.init(cParams);
        end

        function dt = compute(obj,TD)
            x       = obj.designVariable;
            ls2     = x.levelSets{1,2};
            ls3     = x.levelSets{1,3};
            tfiFun  = x.obtainDomainFunction();
            tfi  = obj.splitCellIntoValues(tfiFun);
            chi2 = obj.computeExactCharacteristicFunctionLevelSet(ls2); %- Mixed formulation method
            chi3 = obj.computeExactCharacteristicFunctionLevelSet(ls3); %- Mixed formulation method

            dt(1,:) = - tfi(1,:).*TD{1,end} - tfi(2,:).*TD{2,end} - tfi(3,:).*TD{3,end} ...
                + tfi(4,:).*( (1-chi2).*TD{end,1} + (1-chi3).*chi2.*TD{end,2} + chi2.*chi3.*TD{end,3} );

            dt(2,:) = - tfi(2,:).*TD{2,1} - tfi(3,:).*TD{3,1} + tfi(1,:).*( (1-chi3).*TD{1,2} + chi3.*TD{1,3} );

            dt(3,:) = tfi(2,:).*TD{2,3} - tfi(3,:).*TD{3,2};
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh           = cParams.mesh;
            obj.designVariable = cParams.designVariable;
        end

        function f = computeUnitAnalyticalFunction(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            f         = AnalyticalFunction(s);
        end

        function uFun = computeUnitUnfittedFunction(obj,chiLS)
            s.uMesh = chiLS;
            s.fun   = obj.computeUnitAnalyticalFunction();
            uFun    = UnfittedFunction(s);
        end

        function dV = computeMaxElementVolume(obj)
            q  = Quadrature.create(obj.mesh,2);
            dV = obj.mesh.computeDvolume(q);
            dV = max(sum(dV,1));
        end

        function chi = computeExactCharacteristicFunctionLevelSet(obj,ls)
            chiLS = ls.getUnfittedMesh();
            int   = obj.computeRHSCharFunIntegrator(chiLS);
            uFun  = obj.computeUnitUnfittedFunction(chiLS);
            test  = LagrangianFunction.create(obj.mesh,1,'P0');
            chi   = int.compute(uFun,test);
            dV    = obj.computeMaxElementVolume();
            chi   = chi'/dV;
        end
    end

    methods (Static, Access = private)
        function chiVal = splitCellIntoValues(chi)
            chiVal     = zeros(length(chi),length(chi{1}.fValues));
            for i = 1:length(chi)
                chiVal(i,:) = chi{i}.fValues;
            end
        end

        function int = computeRHSCharFunIntegrator(chiLS)
            s.mesh     = chiLS;
            s.type     = 'Unfitted';
            s.quadType = 2;
            int        = RHSintegrator.create(s);
        end
    end
end