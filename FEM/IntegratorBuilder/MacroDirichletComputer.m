classdef MacroDirichletComputer < handle

    properties (Access = private)
        dirDOFs
        sizeK
    end

    methods (Access = public)
        function obj = MacroDirichletComputer(cParams)
            obj.init(cParams);
        end

        function [CtDir, sizeDir] = computeDirCond(obj)
            [CtDir, sizeDir] = obj.buildDirichletLHS();
        end

    end

    methods (Access = private)
        function init(obj, cParams)
            obj.dirDOFs = cParams.dirDOFs;
            obj.sizeK   = cParams.ndofs;
        end

        function [CtDir, sizeDir] = buildDirichletLHS(obj)
            sizeDir = size(obj.dirDOFs, 1);
            CtDir = zeros(sizeDir, obj.sizeK);
            for i = 1:sizeDir
                CtDir(i,obj.dirDOFs(i)) = 1; 
            end
        end
    end
end