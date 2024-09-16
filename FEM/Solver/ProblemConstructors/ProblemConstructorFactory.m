classdef ProblemConstructorFactory < handle
    
  methods (Access = public, Static)

        function obj = create(cParams)
            switch cParams.type
                case 'Monolithic'
                    switch cParams.mode
                        case 'Disp'
                            obj = solverMonolithicDisp(cParams);
                        case 'Fluc'
                            obj = solverMonolithicFluc(cParams);
                    end
                case 'Reduced'
                    switch cParams.mode
                        case 'Disp'
                            obj = solverReducedDisp(cParams);
                        case 'Fluc'
                            obj = solverReducedFluc(cParams);
                    end
            end
        end

    end

end