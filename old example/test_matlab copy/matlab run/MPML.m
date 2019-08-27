function [mpml_dx,mpml_dz, mpml_dxx, mpml_dzz, mpml_dxx_pzx, mpml_dzz_pxz] = MPML(Vp_max, Node_num ,node_x, node_z, h)


%   z
%%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    %                            %    %
% 14 %             4              % 34 %
%    %                            %    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    %                            %    %
%    %                            %    %
%    %                            %    %
%  1 %                            % 3  %
%    %                            %    %
%    %                            %    %
%    %                            %    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    %                            %    %
% 12 %              2             % 23 %
%    %                            %    %                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% x
%
%%%%%%%%  Set M-PML  %%%%%%%%

delta_Nx = 20;
delta_Nz = 20;
Rcoef    = 0.0001;
USE_MPML_XMIN = true;
USE_MPML_XMAX = false;
USE_MPML_ZMIN = true;
USE_MPML_ZMAX = true;
xmin = min(node_x);
zmin = min(node_z);
xmax = max(node_x);
zmax = max(node_z);

mpml_x_thick = h * delta_Nx;
mpml_z_thick = h * delta_Nz;
xmin_mpml = xmin + mpml_x_thick;
xmax_mpml = xmax - mpml_x_thick;
zmin_mpml = zmin + mpml_z_thick;
zmax_mpml = zmax - mpml_z_thick;

d0_x = 3.0 * Vp_max * log10(1.0/Rcoef) / (2.0 * mpml_x_thick);
d0_z = 3.0 * Vp_max * log10(1.0/Rcoef) / (2.0 * mpml_z_thick);
pxz  = 0.15;
pzx  = 0.15;

mpml_dx = sparse(Node_num,1);
mpml_dz = sparse(Node_num,1);
mpml_dxx = sparse(Node_num,1);
mpml_dzz = sparse(Node_num,1);
mpml_dxx_pzx = sparse(Node_num,1);
mpml_dzz_pxz = sparse(Node_num,1);


for i = 1 : Node_num
    x = node_x(i);
    z = node_z(i);
    
    %%%% area 1 %%%%
    if     ( (x <= xmin_mpml) && (z >= zmin_mpml) && (z <= zmax_mpml) && USE_MPML_XMIN )
        dx      =   d0_x * ( (xmin_mpml-x) / mpml_x_thick )^2 ;
        dz      =   pzx  * dx;
        dxx     = - d0_x * 2 * (xmin_mpml-x) / mpml_x_thick ^2 ;
        dzz     =   0 ;
        dxx_pzx =   pzx  * dxx ;
        dzz_pxz =   pxz  * dzz ;
        
        %%%% area 3 %%%%
    elseif ( (x >= xmax_mpml) && (z >= zmin_mpml) && (z <= zmax_mpml) && USE_MPML_XMAX )
        dx      =   d0_x * ( (x-xmax_mpml) / mpml_x_thick )^2 ;
        dz      =   pzx  * dx ;
        dxx     =   d0_x * 2 * (x-xmax_mpml) / mpml_x_thick ^2 ;
        dzz     =   0 ;
        dxx_pzx =   pzx * dxx ;
        dzz_pxz =   pxz * dzz ;
        
        %%%% area 2 %%%%
    elseif ( (z <= zmin_mpml) && (x >= xmin_mpml) && (x <= xmax_mpml) && USE_MPML_ZMIN )
        dz      =   d0_z * ( (zmin_mpml-z) / mpml_z_thick )^2 ;
        dx      =   pxz  * dz ;
        dxx     =   0 ;
        dzz     = - d0_z * 2 * (zmin_mpml-z) / mpml_z_thick ^2 ;
        dxx_pzx =   pzx  * dxx ;
        dzz_pxz =   pxz  * dzz ;
        
        %%%% area 4 %%%%
    elseif ( (z >= zmax_mpml) && (x >= xmin_mpml) && (x <= xmax_mpml) && USE_MPML_ZMAX )
        dz      =   d0_z * ( (z-zmax_mpml) / mpml_z_thick )^2 ;
        dx      =   pxz  * dz ;
        dxx     =   0 ;
        dzz     =   d0_z * 2 * (z-zmax_mpml) / mpml_z_thick ^2 ;
        dxx_pzx =   pzx  * dxx ;
        dzz_pxz =   pxz  * dzz ;
        
        %%%% area 12 %%%%
    elseif ( (x <= xmin_mpml) && (z <= zmin_mpml) && (USE_MPML_XMIN || USE_MPML_ZMIN ) )
        dx      =   d0_x * ( (xmin_mpml-x) / mpml_x_thick )^2 + pxz * d0_z * ( (zmin_mpml-z) / mpml_z_thick )^2 ;
        dz      =   d0_z * ( (zmin_mpml-z) / mpml_z_thick )^2 + pzx * d0_x * ( (xmin_mpml-x) / mpml_x_thick )^2 ;
        dxx     = - d0_x * 2 * (xmin_mpml-x) / mpml_x_thick ^2 ;
        dzz     = - d0_z * 2 * (zmin_mpml-z) / mpml_z_thick ^2 ;
        dxx_pzx =   pzx  * dxx ;
        dzz_pxz =   pxz  * dzz ;
        
        %%%% area 14 %%%%
    elseif ( (x <= xmin_mpml) && (z >= zmax_mpml) && (USE_MPML_XMIN || USE_MPML_ZMAX ))
        dx      =   d0_x * ( (xmin_mpml-x) / mpml_x_thick )^2 + pxz * d0_z * ( (z-zmax_mpml) / mpml_z_thick )^2 ;
        dz      =   d0_z * ( (z-zmax_mpml) / mpml_z_thick )^2 + pzx * d0_x * ( (xmin_mpml-x) / mpml_x_thick )^2 ;
        dxx     = - d0_x * 2 * (xmin_mpml-x) / mpml_x_thick ^2 ;
        dzz     =   d0_z * 2 * (z-zmax_mpml) / mpml_z_thick ^2 ;
        dxx_pzx =   pzx  * dxx ;
        dzz_pxz =   pxz  * dzz ;
        
        %%%% area 23 %%%%
    elseif ( (x >= xmax_mpml) && (z <= zmin_mpml) && (USE_MPML_XMAX || USE_MPML_ZMIN ))
        
        %     dx      =   d0_x * ( (x-xmax_mpml) / mpml_x_thick )^2 + pxz * d0_z * ( (zmin_mpml-z) / mpml_z_thick )^2 ;
        %     dxx     =   d0_x * 2 * (x-xmax_mpml) / mpml_x_thick ^2 ;
        %     dz      =   d0_z * ( (zmin_mpml-z) / mpml_z_thick )^2 + pzx * d0_x * ( (x-xmax_mpml) / mpml_x_thick )^2 ;
        
        dz      =   d0_z * ( (zmin_mpml-z) / mpml_z_thick )^2 ;            %  modified here
        dx      =   pxz  * d0_z * ( (zmin_mpml-z) / mpml_z_thick )^2  ;    %  modified here
        dxx     =   0.0;                                                   %  modified here
        
        dzz     = - d0_z * 2 * (zmin_mpml-z) / mpml_z_thick ^2 ;
        dxx_pzx =   pzx  * dxx ;
        dzz_pxz =   pxz  * dzz ;
        
        %%% area 34 %%%%
    elseif ( (x >= xmax_mpml) && (z >= zmax_mpml) && (USE_MPML_XMAX || USE_MPML_ZMAX) )
        %     dx      =   d0_x * ( (x-xmax_mpml) / mpml_x_thick )^2 + pxz * d0_z * ( (z-zmax_mpml) / mpml_z_thick )^2 ;
        %     dz      =   d0_z * ( (z-zmax_mpml) / mpml_z_thick )^2 + pzx * d0_x * ( (x-xmax_mpml) / mpml_x_thick )^2 ;
        %     dxx     =   d0_x * 2 * (x-xmax_mpml) / mpml_x_thick ^2 ;
        
        dz      =   d0_z * ( (z-zmax_mpml) / mpml_z_thick )^2  ;          %  modified here
        dx      =   pxz  * d0_z * ( (z-zmax_mpml) / mpml_z_thick )^2 ;    %  modified here
        dxx     =   0.0;                                                  %  modified here
        dzz     =   d0_z * 2 * (z-zmax_mpml) / mpml_z_thick ^2 ;
        dxx_pzx =   pzx  * dxx ;
        dzz_pxz =   pxz  * dzz ;
        
    else
        dx      =   0 ;
        dz      =   0 ;
        dxx     =   0 ;
        dzz     =   0 ;
        dxx_pzx =   pzx  * dxx ;
        dzz_pxz =   pxz  * dzz ;
    end
    
    mpml_dx(i) = dx;
    mpml_dz(i) = dz;
    mpml_dxx(i) = dxx;
    mpml_dzz(i) = dzz;
    mpml_dxx_pzx(i) = dxx_pzx;
    mpml_dzz_pxz(i) = dzz_pxz;

    
end
end
