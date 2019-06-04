classdef VideoMakerFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            
            fileName = cParams.fileName;
            type     = cParams.type;
            dim      = cParams.dim;
            
            if contains(fileName,'SYM','IgnoreCase',true) && contains(fileName,'cantilever','IgnoreCase',true)
                switch dim
                    case '2D'
                        switch type
                            case 'LevelSet'
                                obj = VideoMaker_TopOpt_levelSetBridge();
                            case 'Density'
                                obj = VideoMaker_TopOpt_densityBridge();
                        end
                    case '3D'
                        switch type
                            case 'LevelSet'
                                obj = VideoMaker_TopOpt_levelSet3DYmirror();
                            case 'Density'
                                obj = VideoMaker_TopOpt_density3DYmirror();
                        end
                end
            elseif contains(fileName,'SYM','IgnoreCase',true) && (contains(fileName,'chair','IgnoreCase',true) || contains(fileName,'throne','IgnoreCase',true))
                switch dim
                    case '2D'
                        switch type
                            case 'LevelSet'
                                obj = VideoMaker_TopOpt_levelSetBridge();
                            case 'Density'
                                obj = VideoMaker_TopOpt_densityBridge();
                        end
                    case '3D'
                        switch type
                            case 'LevelSet'
                                obj = VideoMaker_TopOpt_levelSet3DXmirror();
                            case 'Density'
                                obj = VideoMaker_TopOpt_density3DXmirror();
                        end
                end
            elseif contains(fileName,'SYM','IgnoreCase',true)
                switch dim
                    case '2D'
                        switch type
                            case 'LevelSet'
                                obj = VideoMaker_TopOpt_levelSetBridge();
                            case 'Density'
                                obj = VideoMaker_TopOpt_densityBridge();
                        end
                    case '3D'
                        switch type
                            case 'LevelSet'
                                obj = VideoMaker_TopOpt_levelSet3DBridge();
                            case 'Density'
                                obj = VideoMaker_TopOpt_density3DBridge();
                        end
                end
            else
                switch dim
                    case '2D'
                        switch type
                            case 'LevelSet'
                                cParams.tclTemplateName = 'Make_Video_characteristic';                                
                                obj = VideoMaker(cParams);
                            case 'Density'
                                cParams.tclTemplateName = 'Make_Video_density';                                
                                obj = VideoMaker(cParams);
                        end
                    case '3D'
                        switch type
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