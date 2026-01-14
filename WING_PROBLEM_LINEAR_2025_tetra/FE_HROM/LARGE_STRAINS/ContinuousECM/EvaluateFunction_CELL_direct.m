function [DATA] = EvaluateFunction_CELL_direct(DATA,xFE)


% DATA.Integrand =DefaultField(DATA.Integrand,'NumberOfPartitions',1) ;
%  DATA.Integrand.EVALUATE_GRADIENT = 0 ; 
% if DATA.Integrand.NumberOfPartitions > 1
%     DATA.Integrand.EVALUATE_GRADIENT = 0 ;
%     NameParam = DATA.Integrand.NameParameterMatrix ;
%     PARAM = DATA.Integrand.(NameParam) ;
%     if iscell(PARAM)
%         A = cell(size(PARAM)) ;
%         A = A(:)' ;
%         for iparam = 1:length(A)
%             DATA.Integrand.(NameParam) = PARAM{iparam} ;
%             A{iparam} = feval(DATA.Integrand.NameFunctionGenerate,DATA.xLIM,xFE,DATA.Integrand) ;
%         end
%         PARAM = cell2mat(PARAM) ;
%         DATA.Integrand.(NameParam) = PARAM ;
%     else
%         A = feval(DATA.Integrand.NameFunctionGenerate,DATA.xLIM,xFE,DATA.Integrand) ;
%     end
%     
% else
%DATA.Integrand.EVALUATE_GRADIENT = 0 ; 
  %   A = feval(DATA.Integrand.NameFunctionGenerate,DATA.xLIM,xFE,DATA.Integrand) ;
     
     nsnap = size(A,2) ; 
     nrows = floor(nsnap/DATA.Integrand.NumberOfPartitions) ;
    ntimes = floor(nsnap/nrows);
    res =  nsnap-ntimes*nrows;
    if res ~=0
        Partition = [repmat([nrows],ntimes,1);res ] ;
    else
        Partition = repmat([nrows],ntimes,1);
    end
     
     
    % PART= DATA.Integrand.Partition ;
%      ROWS = size(A,1)*ones(1,length(PART)) ; 
%      COLS = PART ; 
 %    A = mat2cell(A',Partition) ; 
  %   A =cellfun(@transpose,A,'UniformOutput',false)' ; 
 ; %end
