function DATAOUT = AverageStressStrain_1st(DATA,stressGLOv,wST,nstrain,ngaus,nelem,DATAOUT,strainGLOv)

%dbstop('4')
if nargin == 0
    load('tmp.mat')
    DATA.NOVOIDS  = 0 ;
end

DATA  = DefaultField(DATA,'VOL_RVE',[]) ;

if  DATA.BCBformulation==1
    % Average stress
    stressGLOv = bsxfun(@times,stressGLOv,wST) ;
    stressGLOv =reshape(stressGLOv,nstrain,[]) ;
    stressAVG = sum(stressGLOv,2) ;
    % Computing volume
    
    % changed 18-JAN-2021
    if isempty(DATA.VOL_RVE)
        ROWS = 1:nstrain:nstrain*ngaus*nelem ;
        wSTs = wST(ROWS);
         VOLUME = sum(wSTs) ;
    else
         VOLUME = DATA.VOL_RVE ;
    end
%     if DATA.NOVOIDS == 1 | isempty(DATA.VOL_RVE)
%         ROWS = 1:nstrain:nstrain*ngaus*nelem ;
%         wSTs = wST(ROWS);
%         VOLUME = sum(wSTs) ;
%     else
%         VOLUME = DATA.VOL_RVE ;
%     end

    
    
    stressAVG = stressAVG/VOLUME ;
    
    DATAOUT.stressAVG = stressAVG   ;
    
    
    strainGLOv = bsxfun(@times,strainGLOv,wST) ;
    strainGLOv =reshape(strainGLOv,nstrain,[]) ;
    strainAVG = sum(strainGLOv,2) ;
    strainAVG = strainAVG/VOLUME ;
    
    disp('------------------------------------------------')
    disp('Difference between input and output average strains')
    try
        diffNORM = (norm(strainAVG-DATA.strainINP)/norm(DATA.strainINP)*100) ;
        disp(['DIFF_e = ',num2str(diffNORM),'%']) ;
        disp('--------------------------------------------')
    catch
    end
    
else
    error('Option not implemented')
end