clc
clear all

load('data.mat')
niter = 1000 ;
SOLVERSM = 2;   %% --> Pre. conj. gradient method 

%%%

%The reverse Cuthill-McKee ordering is obtained with p = symrcm(B);
%R = B(p,p);

if isempty(Gb)
    Kll =  K(DOFl,DOFl) ;
    Fl = F(DOFl) ;
    KlrU = K(DOFl,DOFr)*dR ;
else
    Gast = [Gb,sparse(size(Gb,1),length(DOFf))] ;
    Kll = K(DOFl,DOFl) + (Gast'*K(DOFr,DOFl) + K(DOFl,DOFr)*Gast) + Gast'*(K(DOFr,DOFr)*Gast);
    Fl = [Gb'*F(DOFr) + F(DOFm); F(DOFf)] ;
    KlrU = (Gast'*K(DOFr,DOFr)+K(DOFl,DOFr))*dR ;
end
%
% disp('SOLVING  dL = inv(Kll)* FL')
% dL =Kll\(Fl-KlrU) ;
% aa = toc(aa) ;
% disp(['TIME =',num2str(aa),' s']);
%
% %
aa = tic;

% p = symrcm(Kll) ;
% Kll = Kll(p,p);
% Fl = Fl(p) ;
% KlrU = KlrU(p) ;
if SOLVERSM == 1
    dL =Kll\(Fl-KlrU) ;
elseif SOLVERSM ==2
    dL=    pcg(Kll,Fl-KlrU,[],niter) ;
end
aa = toc(aa) ;
disp(['TIME =',num2str(aa),' s']);