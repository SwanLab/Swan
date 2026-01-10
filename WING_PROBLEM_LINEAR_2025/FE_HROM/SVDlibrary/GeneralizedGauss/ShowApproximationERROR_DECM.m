function  EE =  ShowApproximationERROR_DECM(Xf,z,w,ExactIntegral)


if ~iscell(Xf)
    INTapprox_DECM = Xf(z,:)'*w ; % High-fidelity integral of all snapshots
else
    INTapprox = cell(size(Xf)) ;
    for imat = 1:length(Xf)
        if  ischar(Xf{1})
            load(Xf{imat},'Aloc') ;
            INTapprox{imat} =Aloc(z,:)'*w ;
        else
            INTapprox{imat} =  (Xf{imat}(z,:))'*w ;
        end
    end
    
    INTapprox_DECM = cell2mat(INTapprox') ;
end


disp('---------------------------------------------')
%disp('TOTAL ERROR (SVD of Xf and integration)')

disp('ERROR DECM integration (%)')
EE = norm(ExactIntegral-INTapprox_DECM)/norm(ExactIntegral)*100;
disp([num2str(EE)])
