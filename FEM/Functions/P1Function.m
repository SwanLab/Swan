classdef P1Function < FeFunction

    properties (Access = public)
    end

    properties (Access = private)
        interpolation
    end

    properties (Access = private)
        type
        connec
    end

    methods (Access = public)

        function obj = P1Function(cParams)
            obj.init(cParams);
            obj.createInterpolation();
        end

        function dF = computeDiscontinuousField(obj)
            % Goal: use this function
            fRep = obj.repeatFunctionAtNodes();
            s.fValues = fRep;
            s.connec = obj.connec;
            s.type   = obj.type;
            dF = P1DiscontinuousFunction(s);
        end

        function plot(obj, m) % 2D domains only
            p1DiscFun = obj.computeDiscontinuousField();
            p1DiscFun.plot(m);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.connec = cParams.connec;
            obj.type   = cParams.type;
            obj.fValues = cParams.fValues;
        end

        function createInterpolation(obj)
            m.type = obj.type;
            obj.interpolation = Interpolation.create(m,'LINEAR');
        end

        function fRep = repeatFunctionAtNodes(obj)
           f = obj.fValues;
           nNode  = size(obj.connec,2);
           nDime  = size(f,2);
           nElem  = size(obj.connec,1);
           fNodeElem = zeros(nDime,nNode,nElem);
           fNods  = transpose(f);
           for inode = 1:nNode
               nodes = obj.connec(:,inode);
               fNode = fNods(:,nodes);
               fNodeElem(:,inode,:) = fNode;
           end
           fRep = fNodeElem;
        end

    end

end