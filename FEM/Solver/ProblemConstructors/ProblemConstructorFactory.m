classdef ProblemConstructorFactory < handle
    
  methods (Access = public, Static)

        function obj = create(cParams)
            switch cParams.problemType
                case 'MONOLITHIC'
                    switch cParams.problemMode
                        case 'DISP'
                            obj = ProblemConstructorMonolithicDisp(cParams);
                        case 'FLUC'
                            obj = ProblemConstructorMonolithicFluc(cParams);
                    end
                case 'REDUCED'
                    switch cParams.problemMode
                        case 'DISP'
                            obj = ProblemConstructorReducedDisp(cParams);
                        case 'FLUC'
                            obj = ProblemConstructorReducedFluc(cParams);
                    end
            end
        end

    end

end