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
            s.m              = obj.mesh;
    
            charfun          = CharacteristicFunctionComputer(s); 
            [~,obj.tfi]         = charfun.computeFiandTfi();
        end

        function J = computeFunction(obj)
            surf = obj.area;
            charfunc = obj.tfi;

            V = surf*charfunc';
            for i=1:size(obj.vTar,2)
                J(i) = (V(i)/obj.vTar(i))-1;
            end
            
        end

        function dJ = computeGradient(obj,x)
            n = 3; % size psi
            nnode = obj.mesh.nnodes;
            ls1 = x.designVariable{1,1}.fun.fValues;
            ls2 = x.designVariable{1,2}.fun.fValues;
            ls3 = x.designVariable{1,3}.fun.fValues;

            dJ = zeros(nnode*n,4); % 4 volumes
            dJ(1:nnode,1) = ls2>0;
            dJ(nnode+1:nnode*2,1) = -ls1<=0;
            dJ(1:nnode,2) = ls2<=0 & ls3>0;
            dJ(nnode+1:nnode*2,2) = ls1<=0 & ls3>0;
            dJ(nnode*2+1:nnode*3,2) = -ls1<=0 & ls2>0;
            dJ(1:nnode,3) = ls2<=0 & ls3<=0;
            dJ(nnode+1:nnode*2,3) = ls1<=0 & ls3<=0;
            dJ(nnode*2+1:nnode*3,3) = ls1<=0 & ls2<=0;
            dJ(1:nnode,4) = -1;

            dJ(:,1) = dJ(:,1)/obj.vTar(1);
            dJ(:,2) = dJ(:,2)/obj.vTar(2);
            dJ(:,3) = dJ(:,3)/obj.vTar(3);
            dJ(:,4) = dJ(:,4)/obj.vTar(4);

            dJ(dJ==0) = 1e-6;
        end
    end

    methods (Access = public)
        function title = getTitleToPlot(obj)
            nVol  = size(obj.vTar,2);
            for i = 1:nVol
                title{i,1} = ['Volume',char(string(i))];
            end
        end
    end
end