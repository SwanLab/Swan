function Dcomp = Dmat_kinematicConstraint(fI,BasisUrb,BasisUdef,SingVal_Udef,BasisINTFall,Mdom)

f = cell2mat(fI') ;
BasisUrb_f =BasisUrb(f,:) ;
BasisUdef_f = BasisUdef(f,:) ;

MODIFY_NORM = 0 ; 
if MODIFY_NORM == 1
    % % RB modes with norm = 1
    for i = 1:size(BasisUrb_f,2)
        BasisUrb_f(:,i) = BasisUrb_f(:,i)/norm(BasisUrb_f(:,i)) ;
    end
    for i = 1:size(BasisUdef_f,2)
        BasisUdef_f(:,i) = BasisUdef_f(:,i)/norm(BasisUdef_f(:,i))*SingVal_Udef(i)/SingVal_Udef(1) ;
    end
end

BasisUdom_f = [BasisUrb_f,BasisUdef_f] ;

% STEP 5. Solving minimization problem
if size(BasisINTFall,2) <= size(BasisUdom_f)
    error('Disable option of kinematical constraints')
else
    % WARNING !!! IN FUTURE VERSIONS, ROTATION MATRICES MUST BE
    % INCORPORATED HERE
    Mchol =   chol(Mdom(f,f)) ;
    MBasisINTFall = Mchol*BasisINTFall ;
    MBasisUdom = Mchol*BasisUdom_f ;
    Dcomp = MBasisINTFall\MBasisUdom ;
    
    %    Dcomp = BasisINTFcand\BasisUdom_f ;
end