classdef QRFactorization < Solver

    methods(Static)

        function [Q, R] = householder_qr(LHS,RHS,mesh,bc)
            [m, n] = size(LHS);
            Q = eye(m); 

            for k = 1:min(m-1, n)
                x = LHS(k:m, k);
                e1 = zeros(length(x), 1);
                e1(1) = 1;
                v = sign(x(1)) * norm(x) * e1 + x;
                v = v / norm(v);
                H = eye(m);
                H(k:m, k:m) = eye(length(x)) - 2 * (v * v');
                LHS = H * LHS;
                Q = Q * H';
                QRFactorization.plotQ(Q,mesh,bc,iter)
            end
            
            R = LHS;
            res = norm(LHS*x - RHS);
            QRFactorization.plotQ(res,mesh,bc,iter)
        end

        function plotQ(Q,mesh,bc,iter)
            xFull = bc.reducedToFullVector(Q);
            s.fValues = reshape(xFull,2,[])';
            s.mesh = mesh;
            s.fValues(:,end+1) = 0;
            s.ndimf = 3;
            xF = P1Function(s);
            %xF.plot();
            xF.print(['Q',num2str(iter)],'Paraview')
            fclose('all');
        end

        function plotSolution(res,mesh,bc,iter)
            xFull = bc.reducedToFullVector(res);
            s.fValues = reshape(xFull,2,[])';
            s.mesh = mesh;
            s.fValues(:,end+1) = 0;
            s.ndimf = 3;
            xF = P1Function(s);
            %xF.plot();
            xF.print(['Res',num2str(iter)],'Paraview')
            fclose('all');
        end

    end

end