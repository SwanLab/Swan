function   DATAOUT =   Strength_Composite(DATAOUT,MATERIAL,DATA) ; 

 % ---------------------------------
 % Computing  ULTIMATE STRENGHT AND STRAINS (stress space and space
 % strains)
 % ---------------------------------   

 switch  MATERIAL.FIBER.CONSTITUTIVE_MODEL
     case 'VONMISES'
         switch  MATERIAL.MATRIX.CONSTITUTIVE_MODEL
             case 'VONMISES'
                 DATAOUT =   Strength_Composite_VonMises(DATAOUT,MATERIAL,DATA) ; 
             otherwise
                 error('Option not implemented')
         end
         
     otherwise
         error('Option not implemented')
 end