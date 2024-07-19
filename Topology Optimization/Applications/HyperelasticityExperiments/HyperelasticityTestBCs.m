classdef HyperelasticityTestBCs < handle
    
    methods (Access = public)

        function obj = HyperelasticityBCs(type, perc)
            switch type
                case 'Traction'
                    obj.createBC_2DTraction();
                case 'Hole'
                    obj.createBC_2DHole();
                case 'HoleDirich'
                    obj.createBC_2DHoleDirich(perc);
                case 'Bending'
                    obj.createBC_2DBending();
                case 'Cube'
                    obj.createBC_3DCube();
            end
        end

    end

    methods (Access = private)
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end

end

