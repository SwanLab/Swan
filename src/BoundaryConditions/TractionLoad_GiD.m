classdef TractionLoad_GiD < BoundaryCondition
    
    properties (Access = private)
        type
        values
        fun
    end
    
    methods (Access = public)
        
        function obj = TractionLoad_GiD(mesh,s,type)
            obj.type = type;
            switch type
                case 'DIRAC'
                    if isfield(s,'domain')
                        nodes   = s.domain(mesh.coord);
                        dir     = repmat(s.direction,length(nodes),1);
                        value   = s.value;
                    else
                        nodes   = s.pointload(:,1);
                        dir     = s.pointload(:,2);
                        value   = s.pointload(:,3);
                    end
                    vals = zeros(length(nodes),mesh.ndim);
                    for i = 1:length(value)
                        vals(nodes(i),dir(i)) = value(i);
                    end
                    obj.values = reshape(vals',[],1);
                case 'FUNCTION'
                    dom     = s.domain;
                    f       = s.fun;
                    neuFun  = AnalyticalFunction.create(dom,f.mesh);
                    obj.fun = f.*neuFun;
            end
            obj.mesh = mesh;
        end

        function rhs = computeRHS(obj,test)
            switch obj.type
                case 'DIRAC'
                    rhs = obj.values;
                case 'FUNCTION'
                    f   = obj.fun;
                    rhs = IntegrateRHS(@(v) DP(f,v),test,obj.mesh,'Boundary');
            end
        end
        
    end
end