function [faceDOFS,PhiRB,PsiRBf,Mdom,MdomFF] =   GeometricVar_1dom(REFMESH,DATAROM) 

if nargin == 0
    load('tmp2.mat')
end

nDOM = length(REFMESH);  % Number of domains
% From /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/MultiECMgen/MultiECM_1Dsmall.tex
faceDOFS = cell(1,nDOM) ;
PhiRB = cell(1,nDOM) ;
PsiRBf = cell(1,nDOM) ;
Mdom = cell(1,nDOM) ;
MdomFF = cell(1,nDOM) ;

for  idom = 1:nDOM
    faceDOFS{idom} = cell(2,1) ; 
    faceDOFS{idom}{1} = DATAROM{idom}.f1 ;  % DOFs face 1
    faceDOFS{idom}{2} = DATAROM{idom}.f2 ;  % DOFs face 2
    f1 = DATAROM{idom}.f1 ;  % DOFs face 1
    f2 = DATAROM{idom}.f2 ;  % DOFs face 2
    f = [f1;f2] ;
    PhiRB{idom} = REFMESH{idom}.BasisUrb ;  % Rigid body modes
   
    % Mass matrix subdomain
    Mdom{idom} =REFMESH{idom}.M ;
    % Mass matrix subdomain interfaces     
    MdomFF{idom} = cell(2,2) ;    
    MdomFF{idom}{1,1} = Mdom{idom}(f1,f1);
    MdomFF{idom}{2,2} = Mdom{idom}(f2,f2);
    nrows = size(MdomFF{idom}{1,1},1) ; ncols = size(MdomFF{idom}{2,2},2) ;
    MdomFF{idom}{1,2} = sparse(nrows,ncols) ;
    MdomFF{idom}{2,1} = sparse(ncols,nrows) ;
    
    % Resultant modes 
    PsiRBf{idom} = cell(2,1) ; 
    PsiRBf{idom}{1} = MdomFF{idom}{1,1}* PhiRB{idom}(f1,:) ; 
    PsiRBf{idom}{2} = MdomFF{idom}{2,2}* PhiRB{idom}(f2,:) ; 
    
    
    
    
    
end
