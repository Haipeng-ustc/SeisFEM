function [rho, c11, c13, c33, c44] = VTI_Media ()
%%%%%%%% Elastic VTI Media  %%%%%%%%

% isotropic medium 
rho = 2000.0 ;
scale_aniso = 0.45e9;
c11 = 18.0 * scale_aniso;
c13 = 2.0 * scale_aniso;
c33 = 18.0 * scale_aniso;
c44 = 8.0 * scale_aniso;

% Biotite
% rho = 3050.0 ;
% scale_aniso = 1e9;
% c11 = 27.40 * scale_aniso;
% c13 = 10.53 * scale_aniso;
% c33 = 50.13 * scale_aniso;
% c44 = 5.48  * scale_aniso;

% shale
% rho = 2420.0 ;
% scale_aniso = 1e9;
% c11 = 16.93 * scale_aniso;
% c13 = 14.68 * scale_aniso;
% c33 = 27.60 * scale_aniso;
% c44 = 5.37 * scale_aniso;

% Muscovite 
% rho = 2280.0 ;
% scale_aniso = 1e9;
% c11 = 33.47 * scale_aniso;
% c13 = 11.75 * scale_aniso;
% c33 = 44.54 * scale_aniso;
% c44 = 9.97  * scale_aniso;
% 
% Calcite 
% rho = 2710.0 ;
% scale_aniso = 1e9;
% c11 = 52.48 * scale_aniso;
% c13 = 49.15 * scale_aniso;
% c33 = 77.10 * scale_aniso;
% c44 = 30.47  * scale_aniso;
% 
% Apatite
% rho = 3200.0 ;
% scale_aniso = 1e9;
% c11 = 10.82  * scale_aniso;
% c13 = 59.12  * scale_aniso;
% c33 = 128.63 * scale_aniso;
% c44 = 61.64  * scale_aniso;
% 
% Zinc
% rho = 7100.0 ;
% scale_aniso = 1e9;
% c11 = 165.0  * scale_aniso;
% c13 = 50.0  * scale_aniso;
% c33 = 62.0 * scale_aniso;
% c44 = 39.60  * scale_aniso;



end