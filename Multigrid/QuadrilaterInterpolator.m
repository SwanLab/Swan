classdef QuadrilaterInterpolator < Interpolator

    properties (Access = public)
        I
    end

    properties (Access = private)
        mesh
    end

    methods (Access = public)

        function obj = QuadrilaterInterpolator(cParams)
            obj.init(cParams);
            obj.createInterpolator(obj.mesh);
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end

        function createInterpolator(obj,mesh)
            p = mesh.coord;
            t = mesh.connec;

            n = size(p,1);
            q = size(t,1);
            T = sparse(speye(n,n));
            tnew = []; j = 1;
            p_ori = p;
            for i = 1:q % this will add all the midpoints into p
                tcurr = t(i,:);
                pmid = [
                    (p(tcurr(1),:) + p(tcurr(2),:)) / 2;
                    (p(tcurr(2),:) + p(tcurr(3),:)) / 2;
                    (p(tcurr(3),:) + p(tcurr(4),:)) / 2;
                    (p(tcurr(4),:) + p(tcurr(1),:)) / 2
                    ];
                p = [p; pmid];
                pcenter = (p(5,:) + p(7,:)) / 2;
                p = [p; pcenter];
            end

            [~,ia] = unique(p,'rows','stable');
            Ia = ia(n+1:end);
            ias = ia(n+1:end) - n ;
            potential_tri = ceil(ias./4);
            d = 1;
            midpt_curr = [];

            for i = potential_tri' % now need to loop thru ia and find the triangle that
                %corresponds to this midpoint
                tcurr = t(i,:);
                midpt_curr(1,:) = p(Ia(d),:);

                pmid = [
                    (p(tcurr(1),:) + p(tcurr(2),:)) / 2;
                    (p(tcurr(2),:) + p(tcurr(3),:)) / 2;
                    (p(tcurr(3),:) + p(tcurr(4),:)) / 2;
                    (p(tcurr(4),:) + p(tcurr(1),:)) / 2
                    ];

                if midpt_curr(1,:) == pmid(1,:)
                    T(n + 1, [tcurr(1),tcurr(2)]) = 1/2;
                elseif midpt_curr(1,:) == pmid(2,:)
                    T(n + 1, [tcurr(2),tcurr(3)]) = 1/2;
                elseif midpt_curr(1,:) == pmid(3,:)
                    T(n + 1, [tcurr(3),tcurr(4)]) = 1/2;
                elseif midpt_curr(1,:) == pmid(4,:)
                    T(n + 1, [tcurr(1),tcurr(4)]) = 1/2;
                else
                    T(n + 1, [tcurr(5),tcurr(7)]) = 1/4;
                end
                n = n + 1;
                d = d + 1;
            end
            obj.I = T;
        end
    end
end