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
            for i=1:size(obj.vTar,2)-1
                J(i) = (V(i)/obj.vTar(i))-1;
            end
            
        end

        function dJ = computeGradient(obj,x)
            n = 3; % size psi
            nnode = obj.mesh.nnodes;

            dJ = zeros(nnode*n,3); % 3 unknown volumes
            dJ(1:nnode,1) = 1;
            dJ(nnode+1:nnode*2,1) = -1;
            dJ(1:nnode,2) = 1;
            dJ(nnode+1:nnode*2,2) = 1;
            dJ(nnode*2+1:nnode*3,2) = -1;
            dJ(1:nnode,3) = 1;
            dJ(nnode+1:nnode*2,3) = 1;
            dJ(nnode*2+1:nnode*3,3) = 1;

            dJ(:,1) = dJ(:,1)/obj.vTar(1);
            dJ(:,2) = dJ(:,2)/obj.vTar(2);
            dJ(:,3) = dJ(:,3)/obj.vTar(3);
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