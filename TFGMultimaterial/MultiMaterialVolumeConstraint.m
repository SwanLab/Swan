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
            tfi       = x.obtainDomainFunction(); 
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
            s.mesh           = obj.mesh;
            s.designVariable = x;
            multGrad         = MultimaterialGradientComputer(s);
            TD               = obj.computeTopologicalDerivatives();
            dt               = multGrad.compute(TD);
            dJval            = pdeprtni(obj.mesh.coord',obj.mesh.connec',dt);
            dJval            = reshape(dJval,[],1);
            dJval(dJval==0)  = 1e-6;
            dJ.fValues       = dJval;
        end

        function dV = computeTopologicalDerivatives(obj)
            k  = obj.matID;
            dV = cell(obj.nMat,obj.nMat);
                for i=1:obj.nMat
                    for j=1:obj.nMat
                        if i==j
                            dV{i,j} = zeros(1,obj.mesh.nelem);
                        elseif i == k
                            dV{i,j} = -1/obj.vTar*ones(1,obj.mesh.nelem);
                        elseif j == k
                            dV{i,j} = 1/obj.vTar*ones(1,obj.mesh.nelem);
                        else
                            dV{i,j} = zeros(1,obj.mesh.nelem);
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