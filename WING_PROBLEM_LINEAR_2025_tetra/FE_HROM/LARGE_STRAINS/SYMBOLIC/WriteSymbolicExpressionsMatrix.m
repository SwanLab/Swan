function WriteSymbolicExpressionsMatrix(NameDiary,nstrain,StrInputVariable_symb,StrInputVariable_vect,...
    NameRowsINDEXES,InputVariableSymbolic,OutPutVariableVector)

if nargin == 0
    NameDiary = 'Ivol_diary.txt' ; 
    nstrain = 3; 
    StrInputVariable_symb = 'Cb' ; 
    StrInputVariable_vect = 'Cb' ; 
    NameRowsINDEXES = 'SROWS' ; 
    InputVariableSymbolic = Ivol ; 
    OutPutVariableVector = 'Ivol'  ; 
end

delete(NameDiary)
diary(NameDiary)

INPUT_equiv = cell(nstrain,2) ;
for i=1:nstrain
    INPUT_equiv{i,1} = [StrInputVariable_symb,num2str(i)];
    INPUT_equiv{i,2} = [StrInputVariable_vect,'(',NameRowsINDEXES,'{',num2str(i),'})'];
end

INPUT_var = InputVariableSymbolic;
OUTPUT = OutPutVariableVector ;


for istrain = 1:nstrain
    
    for  jstrain = 1:nstrain
        
        if   jstrain >=istrain
        
        STRLOC = sym2str(INPUT_var(istrain,jstrain)) ;
        
        for hstrain = 1:nstrain
            STRLOC = strrep(STRLOC,INPUT_equiv{hstrain,1},INPUT_equiv{hstrain,2}) ;
        end
        
        Cstr = [OUTPUT,'(',NameRowsINDEXES,'{',num2str(istrain),'}',',',num2str(jstrain),') = ',STRLOC,'; '] ;
        disp(Cstr)
        
        else
            STRLOC = [OUTPUT,'(',NameRowsINDEXES,'{',num2str(jstrain),'}',',',num2str(istrain),')']; 
            Cstr = [OUTPUT,'(',NameRowsINDEXES,'{',num2str(istrain),'}',',',num2str(jstrain),') = ',STRLOC,'; '] ; 
            disp(Cstr)
        end
        
    end
    
end



diary off
