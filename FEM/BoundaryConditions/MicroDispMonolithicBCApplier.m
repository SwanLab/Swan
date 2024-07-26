classdef MicroDispMonolithicBCApplier < handle

    properties (Access = public)
        mesh
        dirichlet_dofs
        vstrain
        ndofs
    end

    methods (Access = public)
        function obj = MicroDispMonolithicBCApplier(cParams)
            obj.init(cParams);
        end

        function [CtDir, CtPerDir] = getLHSMatrix(obj)
            [LDNode, LUNode, RDNode, RUNode] = obj.getVertices();
            [CtDir, CtPerDir]                = obj.computeDirLHS(LDNode,...
                                                LUNode, RDNode, RUNode);
        end

        function [RHSDir, RHSDirPer] = getRHSVector(obj)
            if obj.vstrain(1) == 1 || obj.vstrain(2) == 1
                RHSDir    = zeros(6, 1);
                RHSDirPer = -ones(2,1);
            else
                RHSDir    = zeros(4, 1);
                RHSDir(3) = 0.5;
                RHSDir(4) = 0.5;
                RHSDirPer = -ones(2,1);
            end
        end

    end

    methods (Access = private)
        function init(obj, cParams)
            obj.mesh           = cParams.mesh;
            obj.dirichlet_dofs = cParams.dirichlet_dofs;
            obj.ndofs          = cParams.ndofs;
            obj.vstrain        = cParams.vstrain;
        end

        function [LDNode, LUNode, RDNode, RUNode] = getVertices(obj)
            dirNodes = obj.dirichlet_dofs(2:2:end)/2;
            dirCoord = obj.mesh.coord(dirNodes, :);
            for i = 1:size(dirCoord, 1)
                if dirCoord(i, 1) == 0 && dirCoord(i, 2) == 0
                    LDNode = dirNodes(i);
                elseif dirCoord(i, 1) == 0 && dirCoord(i, 2) ~= 0
                    LUNode = dirNodes(i);
                elseif dirCoord(i, 1) ~= 0 && dirCoord(i, 2) == 0
                    RDNode = dirNodes(i);
                else 
                    RUNode = dirNodes(i);
                end
            end
        end

        function [CtDir, CtPerDir] = computeDirLHS(obj, LDNode, LUNode, RDNode, RUNode)
            if obj.vstrain(1) == 1
                CtDir               = zeros(6, obj.ndofs);
                CtPerDir            = zeros(2, obj.ndofs);
                CtDir(1,LDNode*2-1) = 1;
                CtDir(2,LDNode*2)   = 1;
                CtDir(3,LUNode*2-1) = 1;
                CtDir(4,LUNode*2)   = 1;
                CtDir(5,RDNode*2)   = 1;
                CtDir(6,RUNode*2)   = 1;
                CtPerDir(1,[LDNode*2-1, RDNode*2-1]) = [1 -1];
                CtPerDir(2,[LUNode*2-1, RUNode*2-1]) = [1 -1];
            elseif obj.vstrain(2) == 1
                CtDir               = zeros(6, obj.ndofs);
                CtPerDir            = zeros(2, obj.ndofs);
                CtDir(1,LDNode*2-1) = 1;
                CtDir(2,LDNode*2)   = 1;
                CtDir(3,RDNode*2-1) = 1;
                CtDir(4,RDNode*2)   = 1;
                CtDir(5,LUNode*2-1) = 1;
                CtDir(6,RUNode*2-1) = 1;
                CtPerDir(1,[LDNode*2, LUNode*2]) = [1 -1];
                CtPerDir(2,[RDNode*2, RUNode*2]) = [1 -1];
            else
                CtDir               = zeros(4, obj.ndofs);
                CtPerDir            = zeros(2, obj.ndofs);
                CtDir(1,LDNode*2-1) = 1;
                CtDir(2,LDNode*2)   = 1;
                CtDir(3,RUNode*2-1) = 1;
                CtDir(4,RUNode*2)   = 1;
                CtPerDir(1,[LDNode*2, RDNode*2, LDNode*2-1, LUNode*2-1]) = [1 -1 1 -1];
                CtPerDir(2,[LUNode*2, RUNode*2, RDNode*2-1, RUNode*2-1]) = [1 -1 1 -1];
            end
        end

        function [RHSDir, RHSDirPer] = computeDirRHS(obj)
            if obj.vstrain(1) == 1 || obj.vstrain(2) == 1
                RHSDir    = zeros(6, 1);
                RHSDirPer = -ones(2,1);
            else
                RHSDir    = zeros(4, 1);
                RHSDir(3) = 0.5;
                RHSDir(4) = 0.5;
                RHSDirPer = -ones(2,1);
            end
        end
    
    end
end