classdef MultiMaterialVolumeConstraint < handle

    properties (Access = private)
        mesh
        volumeTarget
        volume
        nMat
        matID
    end

    methods (Access = public)
        function obj = MultiMaterialVolumeConstraint(cParams)
            obj.init(cParams)
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            chi    = obj.computeCharacteristicFunction(x);
            rho    = obj.convertToDensity(chi); % Un altre possibilitat és fer això al tutorial // VolumeTopDer computer per fusionar
            [V,dV] = obj.volume.computeFunctionAndGradient(rho);
            J      = obj.computeFunction(V);
            dJ     = obj.computeGradient(x,dV);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.volumeTarget = cParams.volumeTarget;
            obj.volume       = VolumeFunctional(cParams);
            obj.nMat         = cParams.nMat;
            obj.matID        = cParams.matID;
        end

        function chi = computeCharacteristicFunction(obj,x)
            tfi       = x.obtainGlobalDomainFunction(); 
            k         = obj.matID;
            chiVals   = tfi{k}.fValues;
            s.fValues = chiVals;
            s.mesh    = obj.mesh;
            s.order   = tfi{1}.order;
            chi       = LagrangianFunction(s);
        end

        function rho = convertToDensity(obj,chi)
            s.fun      = chi;
            s.mesh     = obj.mesh;
            s.type     = 'Density';
            s.plotting = false;
            rho        = DesignVariable.create(s);
        end

        function J = computeFunction(obj,V)
            vTar = obj.volumeTarget;
            J    = V/vTar-1;
        end

        function dJ = computeGradient(obj,x,dV)
            TD = obj.computeTopologicalDerivatives(dV{1});
            dt = ChainRule.compute(x,TD);
            dJ = cell(size(dt));
            for i = 1:length(dt)
                dJ{i} = dt{i}.project('P1');
                dJVal = dJ{i}.fValues;
                dJVal(dJVal>=-1e-6 & dJVal<=1e-6) = 1e-6;
                dJ{i}.setFValues(dJVal);
            end
        end

        function tdV = computeTopologicalDerivatives(obj,dV)
            vTar = obj.volumeTarget;
            k    = obj.matID;
            tdV  = cell(obj.nMat,obj.nMat);
            Z    = LagrangianFunction.create(obj.mesh,1,'P1');
                for i=1:obj.nMat
                    for j=1:obj.nMat
                        if i==j
                            tdV{i,j} = Z;
                        elseif i == k
                            tdV{i,j} = -dV./vTar;
                        elseif j == k
                            tdV{i,j} = dV./vTar;
                        else
                            tdV{i,j} = Z;
                        end
                    end
                end
        end
    end

    methods (Access = public)
        function title = getTitleToPlot(obj)
            title = ['Volume mat',char(string(obj.matID))];
        end
    end
end