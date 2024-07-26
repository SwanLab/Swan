classdef PhaseFieldMonitoring < handle

    methods (Access = public, Static)

        function monitor = initialize(cParams)
            s.shallDisplay = cParams.shallDisplay;
            s.mesh = cParams.mesh;
            s.barLim = [0;1];
            switch cParams.type
                case 'Full'
                    s.maxNColumns = 1;
                    s.titles = [{'Force-displacement'},{'Damage'}]; % {'Cost'},{'Energy'}
                    s.chartTypes = [{'plot'},{'surf'}];
                    monitor = Monitoring(s);
                case 'Reduced'
                    s.maxNColumns = 1;
                    s.titles = {'Damage'};
                    s.chartTypes = {'surf'};
                    monitor = Monitoring(s);
            end
        end
    end
    
end