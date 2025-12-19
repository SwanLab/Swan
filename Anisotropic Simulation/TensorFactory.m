
classdef TensorFactory < handle

    methods (Access = public, Static)

        function tensor = create(cParams)

            switch cParams
                case '0'
                    % fileName=fullfile("Elasticity Micro/HomogTensors/C_0.mat");
                    % load(fileName);
                    % tensor=C;
                    tensor = [ 0.4283   0.0012     0;
                               0.0012   0.0035     0;
                               0          0      0.0012];

                case '90'
                    % fileName=fullfile("Elasticity Micro/HomogTensors/C_90.mat");
                    % load(fileName);
                    % tensor=C;
                    tensor = [ 0.0035   0.0012    0;
                               0.0012   0.4283    0;
                                 0      0      0.0012];

                case '45'
                    % fileName=fullfile("Elasticity Micro/HomogTensors/C_45.mat");
                    % load(fileName);
                    % tensor=C;
                    C_0 = [ 0.4283  0.0012     0;
                            0.0012  0.0035     0;
                              0       0      0.0012];
                    theta=-45;

                    Rot = [0.5*(1+cosd(2*theta))  0.5*(1-cosd(2*theta))  sind(2*theta);
                           0.5*(1-cosd(2*theta))  0.5*(1+cosd(2*theta))  -sind(2*theta);
                           -0.5*sind(2*theta)       0.5*sind(2*theta)    cosd(2*theta)];

                    tensor2 = Rot*C_0*Rot';

                    theta=45;
                    Rot = [0.5*(1+cosd(2*theta))  0.5*(1-cosd(2*theta))  sind(2*theta);
                           0.5*(1-cosd(2*theta))  0.5*(1+cosd(2*theta))  -sind(2*theta);
                           -0.5*sind(2*theta)       0.5*sind(2*theta)    cosd(2*theta)];
   
                    tensor = [ 0.1147    0.1127  0.1117;
                                0.1127   0.1147  0.1117;
                                0.1117   0.1117  0.1127];

                    tensor3=Rot*tensor*Rot';

    
                case '0_90'
                    % fileName=fullfile("Elasticity Micro/HomogTensors/C_0_90.mat");
                    % load(fileName);
                    % tensor=C;
                    tensor = [ 0.2825   0.0354  0;
                                0.0354  0.2825  0;
                                0      0      0.0162];

                case '0_45'
                    % fileName=fullfile("Elasticity Micro/HomogTensors/C_0_45.mat");
                    % load(fileName);
                    % tensor=C;
                    tensor = [ 0.2823   0.0615  0.0606;
                               0.0615   0.0826  0.0673;
                               0.0606   0.0673  0.0728];
                case '-45_45'
                    % fileName=fullfile("Elasticity Micro/HomogTensors/C_45_45.mat");
                    % load(fileName);
                    % tensor=C;
                    tensor = [ 0.1701   0.1351  0;
                               0.1351  0.1701 0;
                               0          0      0.1188];
                    
                case '0_45_45_0'
                    tensor = [ 0.5    0.001     0;
                               0.001  0.0022    0;
                               0        0      0.19];
            end

        end

    end
    
end
