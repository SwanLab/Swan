classdef SectionVariableFactory < handle

    methods (Access = public, Static)
        
        function section = create(cParams)
            switch cParams.designVariable.sectionType
                case 'Circular'
                    section = CircularSection(cParams);
                case 'Quadrilateral'
                    section = QuadrilateralSection(cParams);
            end
        end
        
    end

end