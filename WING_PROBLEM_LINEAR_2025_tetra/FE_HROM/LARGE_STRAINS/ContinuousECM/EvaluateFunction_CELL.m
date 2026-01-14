function [A,DATA] = EvaluateFunction_CELL(DATA,xFE,W)

if nargin == 0
    load('tmp.mat')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATA.Integrand.EVALUATE_GRADIENT = 0 ;
DATA = DefaultField(DATA,'LimitMbytesMatricesSnapshots',100) ;
DATA = DefaultField(DATA,'NameWS_Store_Matrices','Amat_') ;

 DATA.Integrand = DefaultField(DATA.Integrand,'NameParameterMatrix',[]) ; 

NameParam = DATA.Integrand.NameParameterMatrix ;  % This should be specified in the input data file
if ~isempty(NameParam)
PARAM = DATA.Integrand.(NameParam) ;
else
    PARAM = [] ; 
end

if ~isempty(PARAM) && iscell(PARAM)
    disp('Generating snapshot matrix (partitioned version)')
    disp(['Block matrices whose size in Mbytes is above ',num2str(DATA.LimitMbytesMatricesSnapshots),' will be stored in memory'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A = cell(size(PARAM)) ;
    A = A(:)' ;
    DATA.ExactIntegral = cell(size(A')); 
    for iparam = 1:length(A)
        disp(['Block iparam = ',num2str(iparam)])
        DATA.Integrand.(NameParam) = PARAM{iparam} ;
        Ai = feval(DATA.Integrand.NameFunctionGenerate,DATA.xLIM,xFE,DATA.Integrand) ;
         DATA.ExactIntegral{iparam} = Ai'*W ; 
        
        nbytesA = prod(size(Ai))*8e-6;
        if nbytesA > DATA.LimitMbytesMatricesSnapshots
            disp(['Storing in memory...(ncols =',num2str(size(Ai,2)),' size in Mb=',num2str(nbytesA), ' )'])
            NameWSLOC = [DATA.NameWS_Store_Matrices,num2str(iparam),'.mat'] ;
            save(NameWSLOC,'Ai') ;
            A{iparam} =NameWSLOC ;
            disp(['...Done '])
        else
            A{iparam}  = Ai;
            disp(['ncols =',num2str(size(Ai,2)),' size in Mb=',num2str(nbytesA), ' '])
        end
        
        
    end
    PARAM = cell2mat(PARAM) ;
    DATA.Integrand.(NameParam) = PARAM ;
    DATA.ExactIntegral = cell2mat(DATA.ExactIntegral) ;  
else
    A = feval(DATA.Integrand.NameFunctionGenerate,DATA.xLIM,xFE,DATA.Integrand) ;
    DATA.ExactIntegral = A'*W ;
    SizeAMb = prod(size(A))*8e-6 ;
    disp(['Matrix Snapshots ',num2str(size(A,1)),' x ',num2str(size(A,2)),'Size  (MBytes) =',num2str(SizeAMb)]) ;
end


% else
%     A = feval(DATA.Integrand.NameFunctionGenerate,DATA.xLIM,xFE,DATA.Integrand) ;
%     DATA.ExactIntegral = A'*W ;
%     SizeAMb = prod(size(A))*8e-6 ;
%     disp(['Matrix Snapshots ',num2str(size(A,1)),' x ',num2str(size(A,2)),'Size  (MBytes) =',num2str(SizeAMb)]) ;
% end
