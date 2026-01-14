function [qLATENT_LOC,qL_extended] = Encoder1MasterMode(BasisU,nREDcoor,SNAP_cluster,DOFl,DATA_evaluateTAU_and_DER,qLATENT_LOC,iloc,idimLAT)


            
            coeff =  BasisU(:,1:nREDcoor)'*SNAP_cluster.DISP.U(DOFl,:) ;
            coeff =bsxfun(@times,coeff',SNAP_cluster.DISP.S)' ;
            qL  = coeff*SNAP_cluster.DISP.V' ;
            [qL_extended,~] = feval(DATA_evaluateTAU_and_DER.nameFunctionEvaluate,qL,DATA_evaluateTAU_and_DER) ;
            qLATENT_LOC{iloc} = qL_extended(idimLAT,:) ;
       