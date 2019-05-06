classdef ConvergenceVarsDispatcher < handle
    
    methods (Access = public, Static)
        
        function names = dispatchNames(optimizer)
            switch optimizer
                case {'AlternatingPrimalDual','DualNestedInPrimal'}
                    names = {'\Deltacost';'Norm L2';'\kappa'};
                case 'MMA'
                    names = {'kktnorm';'outit'};
                case 'IPOPT'
                    names = {'inf_{du}'};
                otherwise
                    error('Invalid optimizer name');
            end
        end
        
    end
end