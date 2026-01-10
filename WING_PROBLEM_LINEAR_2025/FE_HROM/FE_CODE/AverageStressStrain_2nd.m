function DATAOUT = AverageStressStrain_2nd(DATA,stressGLOv,wST,nstrain,ngausLOC,nelem,DATAOUT,strainGLOv,...
    Nst,TypeElement,nnodeE,ndim,COOR,CN,ngaus)

%dbstop('5')
if nargin == 0
    load('tmp.mat')
end

if  DATA.BCBformulation==1
    % Average stress (by integration)
    
    [stressAVG Qx_avg Qy_avg VOLUME]= StressGavg_integration(DATA,stressGLOv,wST,nstrain,ngaus,nelem,DATAOUT,strainGLOv,...
        Nst,TypeElement,nnodeE,ndim,COOR,CN,ngausLOC)
    
    DATAOUT.stressAVG = stressAVG   ;
    DATAOUT.Qy_avg = Qy_avg ;
    DATAOUT.Qx_avg = Qx_avg ;
    
    
    
    
    
     
        strainGLOv = bsxfun(@times,strainGLOv,wST) ;
        strainGLOv =reshape(strainGLOv,nstrain,[]) ;
        strainAVG = sum(strainGLOv,2) ;
        strainAVG = strainAVG/VOLUME ;
%     
%         disp('------------------------------------------------')
%         disp('Difference between input and output average strains')
%         try
%             diffNORM = (norm(strainAVG-DATA.strainINP)/norm(DATA.strainINP)*100) ;
%             disp(['DIFF_e = ',num2str(diffNORM),'%']) ;
%             disp('--------------------------------------------')
%         catch
%         end
    
else
    error('Option not implemented')
end