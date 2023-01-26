classdef SectionVariableFactory < handle

    methods (Access = public, Static)
        
        function section = create(cParams)
            switch cParams.type
                case 'Circular'
                    section = CircularSection(cParams);
                case 'Quadrilateral'
                    section = QuadrilateralSection(cParams);
            end
        end
        
    end

end