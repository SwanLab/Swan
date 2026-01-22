
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

                case '0_3D'
                    tensor=[0.1346 0.0156 0.0154      0       0       0
                            0.0156 0.0578 0.0072      0       0       0
                            0.0154 0.0072 0.0567      0       0       0
                               0      0      0     0.0127     0       0 
                               0      0      0        0    0.0212     0
                               0      0      0        0       0    0.0213];

                case '90_3D'
                    tensor=[0.0446 0.0123 0.0064      0       0       0
                            0.0123 0.1353 0.0072      0       0       0
                            0.0064 0.0072 0.0445      0       0       0
                               0      0      0     0.0178     0       0 
                               0      0      0        0    0.0116     0
                               0      0      0        0       0    0.0178];

                case '45_3D'
                    tensor=[0.0296 0.0275 0.0056      0          0    0.0259
                            0.0275 0.0296 0.0056      0          0    0.0259   
                            0.0056 0.0056 0.0383      0          0    0.0044
                               0      0      0     0.0082     0.0071     0 
                               0      0      0     0.0071     0.0082     0
                            0.0259 0.0259 0.0044      0          0    0.0271];


                case '0_90_3D'
                    tensor=[0.0913 0.0150 0.0110      0        0         0
                            0.0150 0.0995 0.0106      0        0         0
                            0.0110 0.0106 0.0467      0    -0.0022       0
                               0      0      0     0.0160      0      -0.0013 
                               0      0  -0.0022      0     0.0212       0
                               0      0      0    -0.0013        0    0.0213];

                case '0_45_3D'
                    tensor=[ 0.0854  0.0218   0.0081     0    -0.0016  0.0130;
                             0.0218  0.0466   0.0060  0.0013      0    0.0129;
                             0.0081  0.0060   0.0329  0.0023  -0.0025  0.0013;
                                0    0.0013   0.0023  0.0091   0.0035     0;
                             -0.0016   0     -0.0025  0.0035   0.0100     0;
                             0.0130  0.0129   0.0013     0        0    0.0254 ];

                case '45_45_3D'
                    tensor=[ 0.0332  0.0243  0.0052     0    -0.0013    0.0033;
                            0.0243  0.0326   0.0044     0        0      0.0033;
                            0.0052  0.0044   0.0227     0    -0.0027       0;
                               0       0        0    0.0057   0.0010       0;
                           -0.0013     0    -0.0027  0.0010   0.0060       0;
                            0.0033  0.0033      0       0        0      0.0256 ];

            end

        end

    end
    
end
