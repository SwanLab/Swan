function baselineDerivQuality(paramB,paramRho,C,fun,dfun,tag)
if nargin < 6, tag = 'baseline'; end
pB = paramB(:)';  pR = paramRho(:)';
[BB,RR] = meshgrid(pB,pR);                        % nRho x nB, igual a C(...,iRho,iB)
in = false(size(BB));  in(3:end-2,3:end-2) = true; % exclui bordas

ids   = {[1 1 1 1],[1 1 2 2],[1 1 1 2],[2 2 2 2],[2 2 1 2],[1 2 1 2]};
names = {'C1111','C1122','C1112','C2222','C2212','C1212'};

fprintf('\n===== DERIVATIVE QUALITY: %s =====\n',tag);
fprintf('%-7s %10s %10s %9s %11s %10s\n', ...
    'comp','nRMSE(C)','nRMSE db','sign db','nRMSE drho','sign drho');
for ic = 1:6
    id = ids{ic};
    Cg = squeeze(C(id(1),id(2),id(3),id(4),:,:));
    [dRb,dRr] = gradient(Cg,pB,pR);

    v = fun{id(1),id(2),id(3),id(4)}(BB,RR);
    d = dfun{id(1),id(2),id(3),id(4)}(BB,RR);

    eV  = norm(v(in)-Cg(in))/norm(Cg(in));
    ref = {dRb,dRr};  out = zeros(1,4);
    for q = 1:2
        r = ref{q}(in);  n = d{q}(in);
        relev = abs(r) > 0.05*max(abs(r));
        out(2*q-1) = norm(n-r)/norm(r);
        out(2*q)   = 100*mean(sign(n(relev)) ~= sign(r(relev)));
    end
    fprintf('%-7s %10.3e %10.3e %8.2f%% %11.3e %9.2f%%\n', ...
        names{ic},eV,out(1),out(2),out(3),out(4));
end
fprintf('escala max|C| = %.4e\n',max(abs(C(:))));
end