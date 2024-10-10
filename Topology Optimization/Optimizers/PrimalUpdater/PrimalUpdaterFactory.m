classdef PrimalUpdaterFactory < handle
    
   
    methods (Access = public, Static)
    
        function p = create(cParams)
            switch cParams.primal
                case 'SLERP'
                    ls = cParams.designVariable;
                    switch class(ls)
                        case 'LevelSet'
                            s.mesh = ls.fun.mesh;
                        case 'MultiLevelSet'
                            s.mesh = ls.levelSets{1,1}.fun.mesh;
                    end
                    p = SLERP(s);
                case 'PROJECTED GRADIENT'
                    p = ProjectedGradient(cParams);
                case 'HAMILTON-JACOBI'
                    p = HamiltonJacobi(cParams);
            end
        end
    
    end
    
end