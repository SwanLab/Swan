classdef VideoMaker_TopOpt < VideoMaker_Physical_Problem
    
    
    properties
        
        
    end
    
    
    methods (Access = public)
        function obj = VideoMaker_TopOpt()
            
            
        end
    end
    
    methods (Static)
        
        function obj = Create(optimizer,pdim,case_file)
            if contains(case_file,'SYM','IgnoreCase',true) && contains(case_file,'cantilever','IgnoreCase',true)
                switch pdim
                    case '2D'
                        switch optimizer
                            case {'SLERP', 'HAMILTON-JACOBI'}
                                obj = VideoMaker_TopOpt_levelSetBridge();
                            case 'PROJECTED GRADIENT'
                                obj = VideoMaker_TopOpt_densityBridge();
                            case 'MMA'
                                obj = VideoMaker_TopOpt_densityBridge();
                            case 'IPOPT'
                                obj = VideoMaker_TopOpt_densityBridge();
                        end
                    case '3D'
                        switch optimizer
                            case {'SLERP', 'HAMILTON-JACOBI'}
                                obj = VideoMaker_TopOpt_levelSet3DYmirror();
                            case 'PROJECTED GRADIENT'
                                obj = VideoMaker_TopOpt_density3DYmirror();
                            case 'MMA'
                                obj = VideoMaker_TopOpt_density3DYmirror();
                            case 'IPOPT'
                                obj = VideoMaker_TopOpt_density3DYmirror();
                        end
                end
                
            elseif contains(case_file,'SYM','IgnoreCase',true) && (contains(case_file,'chair','IgnoreCase',true) || contains(case_file,'throne','IgnoreCase',true))
                switch pdim
                    case '2D'
                        switch optimizer
                            case {'SLERP', 'HAMILTON-JACOBI'}
                                obj = VideoMaker_TopOpt_levelSetBridge();
                            case 'PROJECTED GRADIENT'
                                obj = VideoMaker_TopOpt_densityBridge();
                            case 'MMA'
                                obj = VideoMaker_TopOpt_densityBridge();
                            case 'IPOPT'
                                obj = VideoMaker_TopOpt_densityBridge();
                        end
                    case '3D'
                        switch optimizer
                            case {'SLERP', 'HAMILTON-JACOBI'}
                                obj = VideoMaker_TopOpt_levelSet3DXmirror();
                            case 'PROJECTED GRADIENT'
                                obj = VideoMaker_TopOpt_density3DXmirror();
                            case 'MMA'
                                obj = VideoMaker_TopOpt_density3DXmirror();
                            case 'IPOPT'
                                obj = VideoMaker_TopOpt_density3DXmirror();
                        end
                end
            elseif contains(case_file,'SYM','IgnoreCase',true)
                switch pdim
                    case '2D'
                        switch optimizer
                            case {'SLERP', 'HAMILTON-JACOBI'}
                                obj = VideoMaker_TopOpt_levelSetBridge();
                            case 'PROJECTED GRADIENT'
                                obj = VideoMaker_TopOpt_densityBridge();
                            case 'MMA'
                                obj = VideoMaker_TopOpt_densityBridge();
                            case 'IPOPT'
                                obj = VideoMaker_TopOpt_densityBridge();
                        end
                    case '3D'
                        switch optimizer
                            case {'SLERP', 'HAMILTON-JACOBI'}
                                obj = VideoMaker_TopOpt_levelSet3DBridge();
                            case 'PROJECTED GRADIENT'
                                obj = VideoMaker_TopOpt_density3DBridge();
                            case 'MMA'
                                obj = VideoMaker_TopOpt_density3DBridge();
                            case 'IPOPT'
                                obj = VideoMaker_TopOpt_density3DBridge();
                        end
                end
            else
                
                switch pdim
                    case '2D'
                        switch optimizer
                            case {'SLERP', 'HAMILTON-JACOBI'}
                                obj = VideoMaker_TopOpt_levelSet();
                            case 'PROJECTED GRADIENT'
                                obj = VideoMaker_TopOpt_density();
                            case 'MMA'
                                obj = VideoMaker_TopOpt_density();
                            case 'IPOPT'
                                obj = VideoMaker_TopOpt_density();
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
                        end
                end
            end
        end
    end
    
end