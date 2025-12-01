classdef TensorFactory < handle

    methods (Access = public, Static)

        function tensor = create(cParams)

            switch cParams
                case '0'
                    tensor = [ 0.4283   0.0012  0;
                               0.0012  0.0035 0;
                               0      0      0.0048];

                case '90'
                    tensor = [ 0.0035   0.0012  0;
                            0.0012  0.4283 0;
                            0      0      0.0048];

                case '45'
                    tensor = [ 0.1147   0.1127  0.2234;
                               0.1127  0.1147 0.2234;
                               0.2234      0.2234   0.4507];
    
                case '0_90'
                    tensor = [ 0.2825   0.0354  0;
                                0.0354  0.2825  0;
                                0      0      0.0648];

                case '0_45'
                    tensor = [ 0.2823   0.0615  0.1212;
                               0.0615  0.0826 0.1346;
                               0.1212      0.1346      0.2913];
                case '-45_45'
                    tensor = [ 0.1701   0.1351  0;
                               0.1351  0.1701 0;
                               0      0      0.4752];
                case '0_45_45_0'
                    tensor = [ 0.5   0.001  0;
                               0.001  0.0022 0;
                               0      0      0.19];
            end

        end

    end
    
end
