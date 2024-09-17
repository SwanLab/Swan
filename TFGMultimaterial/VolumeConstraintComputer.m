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
            s.designVariable = x;
            s.mesh           = obj.mesh;
    
            charfun          = MultiMaterialCharacteristicFunction(s); 
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
            s.designVariable = x;
            s.mesh           = obj.mesh;

            charfun = MultiMaterialCharacteristicFunction(s);
            [~,tfiFun] = charfun.computeAtNodesAndElements();
            tfi = tfiFun.fValues';
            psi2 = x.levelSets{1,2}.fun.fValues;
            psi3 = x.levelSets{1,3}.fun.fValues;

            t = obj.mesh.connec';
            p = obj.mesh.coord';
            [tXi2,~] = integ_exact(t,p,psi2); chi2 = (1 - tXi2); %- Mixed formulation method
            [tXi3,~] = integ_exact(t,p,psi3); chi3 = (1 - tXi3); %- Mixed formulation method
            %     fi = (pdeintrp(p,t,fi)).'; % interpolation at gauss point - P1 projection method
            %     chi2 = (pdeintrp(p,t,(psi(:,2)<0))).'; %- P1 projection method

            dt = [];
            dt(1,:) = - tfi(1,:).*TD{1,end} - tfi(2,:).*TD{2,end} - tfi(3,:).*TD{3,end} ...
                + tfi(4,:).*( (1-chi2).*TD{end,1} + (1-chi3).*chi2.*TD{end,2} + chi2.*chi3.*TD{end,3} );

            dt(2,:) = - tfi(2,:).*TD{2,1} - tfi(3,:).*TD{3,1} + tfi(1,:).*( (1-chi3).*TD{1,2} + chi3.*TD{1,3} );

            dt(3,:) = tfi(2,:).*TD{2,3} - tfi(3,:).*TD{3,2};
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