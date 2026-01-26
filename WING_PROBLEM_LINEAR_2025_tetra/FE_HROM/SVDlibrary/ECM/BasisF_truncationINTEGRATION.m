function   [BasisF] = BasisF_truncationINTEGRATION(intEXACT,U,S,V,Wdom,DATAIN) 

if nargin == 0
    load('tmp.mat')
end

disp('---------------------------------------------------------------')
disp(['Truncation criterion based on integration error '])
disp('---------------------------------------------------------------')

%DATAIN = DefaultField(DATAIN,'TOL_LOC_InternalForces',1e-6) ; 
TOL = DATAIN.TOL_LOC_InternalForces ; % Relative tolerance integration error  

US  = bsxfun(@times,U',S)' ; 

intAPPROX = zeros(1,size(V,1)) ; 
nintEXACT = norm(intEXACT) ; 
errorGLO = zeros(length(S),1) ; 
errorFRO = zeros(length(S),1) ; 
nerrorFRO  = sqrt(sum(S.^2)); 
eREL = 1e20 ; 

imode = 1; 

while eREL > TOL 
    intAPPROX = intAPPROX + sqrt(Wdom)'*(US(:,imode)*V(:,imode)') ; 
    eREL =  norm(intAPPROX-intEXACT)/nintEXACT ; 
    eRELfro = sqrt(sum(S(imode+1:end).^2))/nerrorFRO ; 
    errorFRO(imode) = eRELfro ; 
    errorGLO(imode) = eREL ; 
    disp(['Mode = ',num2str(imode),' (of ',num2str(length(S))...
        ,'): Int. error (%)= ',num2str(eREL*100),...
        '; Frob. error % = ',num2str(eRELfro*100)]); 

    imode = imode + 1; 
end

BasisF = U(:,1:imode-1) ; 

figure(679)
hold on 
xlabel('Mode')
ylabel('log( error (rel))')
h1 = plot(log(errorGLO(1:imode-1)),'Color',rand(1,3)) 
h2 = plot(log(errorFRO(1:imode-1)),'Color',rand(1,3)) 
legend([h1,h2],{'Integ. error','Frob. error'})