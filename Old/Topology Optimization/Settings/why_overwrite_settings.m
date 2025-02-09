clc
clear all

default=SettingsNumericalHomogenizer();
fprintf('Default type: %s\n',default.levelSetCreatorParams.levelSetType)
clear
a=Settings_TestNumericalHomogenizerCustom();
fprintf('Custom type: %s\n',a.levelSetCreatorParams.levelSetType)
b=SettingsNumericalHomogenizer();
fprintf('New default instance after creating custom type: %s\n',b.levelSetCreatorParams.levelSetType)
