classdef TensorFactory < handle

    methods (Access = public, Static)

        function tensor = create(cParams)

            switch cParams
                case '0'
                    tensor = [ 0.4283   0.0012  0;
                               0.00012  0.0035 0;
                               0      0      0.0048];

                case '90'
                    tensor = [ 0.0022   0.001  0;
                                0.001    0.5   0;
                                0        0     0.19];

                case '45'
                    C_0 = [ 0.5   0.001  0;
                            0.001  0.0022 0;
                            0      0      0.19];
                    Rot = [cos(45) -sin(45) 0;
                            sin(45) cos(45) 0;
                            0         0     1];

                    tensor = Rot'*C_0*Rot;
                case '0_90'
                    tensor = [ 0.2076   0.0171  0;
                                0.0171  0.2076  0;
                                0      0      0.0175];
                case '0_45'
                    tensor = [ 0.5   0.001  0;
                               0.001  0.0022 0;
                               0      0      0.19];
                case '-45_45'
                    tensor = [ 0.5   0.001  0;
                               0.001  0.0022 0;
                               0      0      0.19];
                case '0_45_45_0'
                    tensor = [ 0.5   0.001  0;
                               0.001  0.0022 0;
                               0      0      0.19];
            end

        end

    end
    
end
