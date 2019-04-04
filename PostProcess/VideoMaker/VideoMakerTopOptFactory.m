classdef VideoMakerTopOptFactory < handle
    
    methods (Access = public)
        
        function obj = VideoMakerTopOptFactory()
            
        end
        
        function obj = create(obj,filename,varType,dim)
            if contains(filename,'SYM','IgnoreCase',true) && contains(filename,'cantilever','IgnoreCase',true)
                switch dim
                    case '2D'
                        switch varType
                            case 'LevelSet'
                                obj = VideoMaker_TopOpt_levelSetBridge();
                            case 'Density'
                                obj = VideoMaker_TopOpt_densityBridge();
                        end
                    case '3D'
                        switch varType
                            case 'LevelSet'
                                obj = VideoMaker_TopOpt_levelSet3DYmirror();
                            case 'Density'
                                obj = VideoMaker_TopOpt_density3DYmirror();
                        end
                end
            elseif contains(filename,'SYM','IgnoreCase',true) && (contains(filename,'chair','IgnoreCase',true) || contains(filename,'throne','IgnoreCase',true))
                switch dim
                    case '2D'
                        switch varType
                            case 'LevelSet'
                                obj = VideoMaker_TopOpt_levelSetBridge();
                            case 'Density'
                                obj = VideoMaker_TopOpt_densityBridge();
                        end
                    case '3D'
                        switch varType
                            case 'LevelSet'
                                obj = VideoMaker_TopOpt_levelSet3DXmirror();
                            case 'Density'
                                obj = VideoMaker_TopOpt_density3DXmirror();
                        end
                end
            elseif contains(filename,'SYM','IgnoreCase',true)
                switch dim
                    case '2D'
                        switch varType
                            case 'LevelSet'
                                obj = VideoMaker_TopOpt_levelSetBridge();
                            case 'Density'
                                obj = VideoMaker_TopOpt_densityBridge();
                        end
                    case '3D'
                        switch varType
                            case 'LevelSet'
                                obj = VideoMaker_TopOpt_levelSet3DBridge();
                            case 'Density'
                                obj = VideoMaker_TopOpt_density3DBridge();
                        end
                end
            else
                switch dim
                    case '2D'
                        switch varType
                            case 'LevelSet'
                                obj = VideoMaker_TopOpt_levelSet();
                            case 'Density'
                                obj = VideoMaker_TopOpt_density();
                        end
                    case '3D'
                        switch varType
                            case 'LevelSet'
                                obj = VideoMaker_TopOpt_levelSet3D();
                            case 'Density'
                                obj = VideoMaker_TopOpt_density3D();
                        end
                end
            end
        end
        
    end
    
end