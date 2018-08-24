classdef VideoMaker_TopOpt < VideoMaker_Physical_Problem
    methods (Static)
        function obj = Create(optimizer,pdim)
            switch pdim
                case '2D'
                    switch optimizer
                        case {'SLERP', 'HAMILTON JACOBI'}
                            obj = VideoMaker_TopOpt_levelSet();
                        case 'PROJECTED GRADIENT'
                            obj = VideoMaker_TopOpt_density();
                        case 'MMA'
                            obj = VideoMaker_TopOpt_density();
                        case 'IPOPT'
                            obj = VideoMaker_TopOpt_density();
                        otherwise
                            error('Invalid optimizer');
                    end
                case '3D'
                    switch optimizer
                        case {'SLERP', 'HAMILTON-JACOBI'}
                            obj = VideoMaker_TopOpt_levelSet3D();
                        case 'PROJECTED GRADIENT'
                            obj = VideoMaker_TopOpt_density3D();
                        case 'MMA'
                            obj = VideoMaker_TopOpt_density3D();
                        case 'IPOPT'
                            obj = VideoMaker_TopOpt_density3D();
                        otherwise
                            error('Invalid optimizer');
                    end
            end
        end
    end
end