classdef FilterFactory < handle

    methods (Access = public, Static)

        function filter = create(cParams)
            switch cParams.filterType
                case 'P1'
                    filter = FilterKernel(cParams);
                case 'PDE'
                    if not(isfield(cParams,'boundaryType'))
                        cParams.boundaryType = 'Neumann';
                    end
                    if not(isfield(cParams,'metric'))
                        cParams.metric       = 'Isotropy';
                    end
                    switch cParams.boundaryType
                        case {'Neumann','Periodic'}
                            switch cParams.metric
                                case 'Isotropy'
                                    cParams.LHSint.domain = @(e,u,v) e^2.*DP(Grad(v),Grad(u)) + DP(v,u);
                                case 'Anisotropy'
                                    A = cParams.A;
                                    cParams.LHSint.domain = @(e,u,v) e^2.*DP(Grad(v),DP(A,Grad(u)',2,1),2,1) + DP(v,u);
                            end
                            cParams.LHSint.boundary = [];
                        case 'Robin'
                            switch cParams.metric
                                case 'Isotropy'
                                    cParams.LHSint.domain = @(e,u,v) e^2.*DP(Grad(v),Grad(u)) + DP(v,u);
                                case 'Anisotropy'
                                    A = cParams.A;
                                    cParams.LHSint.domain = @(e,u,v) e^2.*DP(Grad(v),DP(A,Grad(u)',2,1),2,1) + DP(v,u);
                            end
                            cParams.LHSint.boundary = @(e,u,v) e.*DP(v,u);
                        case 'DirichletProjection'
                            cParams.LHSint.domain = @(e,u,v) DP(v,u);
                            cParams.LHSint.boundary = @(e,u,v) e.*DP(v,u);
                    end
                    filter = FilterPDE(cParams);
                case 'LUMP'
                    filter = FilterLump(cParams);
                case 'FilterAndProject'
                    filter = FilterAndProject(cParams);
                case 'FilterAdjointAndProject'
                    filter = FilterAdjointAndProject(cParams);
                case 'CloseOperator'
                    filter = CloseOperator(cParams);
                case 'CloseAdjointOperator'
                    filter = CloseAdjointOperator(cParams);
                case 'Segment'
                    filter = NonLinearFilterSegment(cParams);
                case 'Droplet'
                    filter = NonLinearFilterDroplet(cParams);
            end
        end

    end

end