classdef VolumeConstraintComputer < handle
    
    properties (Access = public)
        J
        dJ
    end

    properties (Access = private)
        area
        mesh
        vTar
        tfi
    end

    methods (Access = public)

        function obj = VolumeConstraintComputer(cParams)
            obj.init(cParams)
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            obj.computeCharacteristicFunction(x);
            J = obj.computeFunction();
            dJ = obj.computeGradient(x);
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.area  = cParams.area;
            obj.vTar = cParams.volumeTarget;
        end

        function computeCharacteristicFunction(obj,x)
            charfun          = x.obtainDomainFunction(); 
            [~,obj.tfi]         = charfun.computeAtNodesAndElements();
        end

        function J = computeFunction(obj)
            surf = obj.area;
            charfunc = obj.tfi.fValues';

            V = surf*charfunc';
            for i=1:size(obj.vTar,2)-1
                J(i) = (V(i)/obj.vTar(i))-1;
            end
            
        end

        function dJ = computeGradient(obj,x)
            n = 3; % size psi
            for k=1:n
            dC = cell(n+1,n+1);


            for i=1:(n+1)
                for j=1:(n+1)
                    if i==j
                        dC{i,j} = zeros(1,obj.mesh.nelem);
                    elseif i == k
                        dC{i,j} = -1/obj.vTar(k)*ones(1,obj.mesh.nelem);
                    elseif j == k
                        dC{i,j} = 1/obj.vTar(k)*ones(1,obj.mesh.nelem);
                    else
                       dC{i,j} = zeros(1,obj.mesh.nelem);
                  %      dC{i,j,k} = 1/obj.vTar(j)-1/obj.vTar(i);
                    end
                end
            end
            dt = obj.smoothingAndChainRuleComputing(dC,x);
            t = obj.mesh.connec';
            p = obj.mesh.coord';
            dJk = pdeprtni(p,t,dt);
            dJk = reshape(dJk,[],1);
            dJ(:,k) = dJk;


            end





            % intent
            % dt1 = dt.*Tfi(1,:);
            % dt2 = dt.*Tfi(2,:);
            % dt3 = dt.*Tfi(3,:);
            % 
            % dJ1 = pdeprtni(p,t,dt1);
            % dJ2 = pdeprtni(p,t,dt2);
            % dJ3 = pdeprtni(p,t,dt3);
            % 
            % dJ1 = reshape(dJ1,[],1);
            % dJ2 = reshape(dJ2,[],1);
            % dJ3 = reshape(dJ3,[],1);

            %dJ = reshape(dJ,[],1);
            %dJ = [dJ1, dJ2, dJ3];
            dJ(dJ==0) = 1e-6;

        end

        function dt = smoothingAndChainRuleComputing(obj,TD,x)
            s.mesh           = obj.mesh;
            s.designVariable = x;
            multGrad         = MultimaterialGradientComputer(s);
            dt               = multGrad.compute(TD);
        end
    end

    methods (Access = public)
        function title = getTitleToPlot(obj)
            nVol  = size(obj.vTar,2);
            for i = 1:nVol-1
                title{i,1} = ['Volume',char(string(i))];
            end
        end
    end
end