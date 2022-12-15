classdef SectionVariableFactory < handle

    methods (Access = public, Static)
        
        function section = create(cParams)
            switch cParams.type
                case 'Circular'
                    section = CircularSection(cParams);
            end
        end
        
    end

end