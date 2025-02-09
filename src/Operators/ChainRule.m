classdef ChainRule < handle

    methods(Static, Access=public)

        function dJDV = compute(x,dJMat)
            switch class(x)
                case {'MultiLevelSet'}
                    chR  = MultimaterialGradientComputer(x);
                    dJDV = chR.compute(dJMat);
                otherwise
                    dJDV = dJMat;
            end
        end
    end
end