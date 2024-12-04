classdef MultiMaterialVolumeConstraint < handle

    properties (Access = private)
        mesh
        vTar
        nMat
        matID
    end

    methods (Access = public)

        function obj = MultiMaterialVolumeConstraint(cParams)
            obj.init(cParams)
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            chi = obj.computeCharacteristicFunction(x);
            J   = obj.computeFunction(chi);
            dJ  = obj.computeGradient(x);
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.vTar  = cParams.volumeTarget;
            obj.nMat  = cParams.nMat;
            obj.matID = cParams.matID;
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

        function J = computeFunction(obj,chi)
            quad     = Quadrature.create(obj.mesh,2);
            dV       = obj.mesh.computeDvolume(quad);
            totalVol = sum(dV(:));
            V        = Integrator.compute(chi,obj.mesh,quad.order);
            V        = V/totalVol;
            J        = (V/obj.vTar)-1;
        end

        function dJ = computeGradient(obj,x)
            multGrad         = MultimaterialGradientComputer(x);
            TD               = obj.computeTopologicalDerivatives();
            dt               = multGrad.compute(TD);
            for i = 1:length(dt)
                dJ{i} = dt{i}.project('P1',obj.mesh);
                dJ{i}.fValues(dJ{i}.fValues>=-1e-6 & dJ{i}.fValues<=1e-6) = 1e-6;
            end
        end

        function dV = computeTopologicalDerivatives(obj)
            k  = obj.matID;
            dV = cell(obj.nMat,obj.nMat);
            Z  = LagrangianFunction.create(obj.mesh,1,'P1');
            I  = LagrangianFunction.create(obj.mesh,1,'P1');
            I.fValues = ones(size(I.fValues));
                for i=1:obj.nMat
                    for j=1:obj.nMat
                        if i==j
                            dV{i,j} = Z;
                        elseif i == k
                            dV{i,j} = -I./obj.vTar;
                        elseif j == k
                            dV{i,j} = I./obj.vTar;
                        else
                            dV{i,j} = Z;
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