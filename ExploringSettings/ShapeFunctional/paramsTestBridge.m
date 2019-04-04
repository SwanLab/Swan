filterParams = SettingsFilter('paramsFilterP1HJ');
domainType = 'MACRO';
materialInterpolationParams = SettingsInterpolation('paramsMatInterpSIMPP3ISO2D');
materialInterpolationParams.constitutiveProperties.nu_plus   = 0.3;
materialInterpolationParams.constitutiveProperties.nu_minus  = 0.3;
filename = 'Bridge_Quadrilateral_Bilinear_Structured';