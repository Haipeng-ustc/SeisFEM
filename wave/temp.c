   /********************************************************************************************************************************************
        Equation u1:
         mass * U1tt_new = - c11 * dphidx * dphidx * U_now - 2.0 * rho * mpml_dx * phi * phi * U1t_now - rho * mpml_dx * mpml_dx * phi * phi * U1_now 
                           + phi * phi * Lx1_now + phi * phi * Lx2_now + Source_x
        *********************************************************************************************************************************************/
		for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; y3[i] = 0.0; y4[i] = 0.0; y5[i] = 0.0;
        }
        csr_matvec(csr_p_size, stif1_csr_p, stif1_csr_j, stif1_csr_x,   U_now, y1); // dphidx * dphidx *  U_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, U1t_now, y2); // phi    * phi    * U1t_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x,  U1_now, y3); // phi    * phi    * U1_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Lx1_now, y4); // phi    * phi    * Lx1_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Lx1_now, y5); // phi    * phi    * Lx2_now
        for (i = 0; i < node_num; i++)
        {
            rhs1[i]  =  - c[0][i] * y1[i] - 2.0 * rho[i] * mpml_dx[i] * y2[i] -  rho[i] * mpml_dx[i] * mpml_dx[i] * y3[i] 
                        + y4[i] + y5[i] +  (i == source_node) * point_source * sin( Angle_force * pi / 180.0 );
        }
       
       /********************************************************************************************************************************************
        Equation u2:
         mass * U2tt_new = - c13 * dphidx * dphidy * W_now - c44 * dphidy * dphidx * W_now - rho * mpml_dx * phi * phi * U2t_now 
                           - rho * mpml_dy * phi * phi * U2t_now - rho * mpml_dx * mpml_dy * phi * phi * U2_now
        *********************************************************************************************************************************************/
        for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; y3[i] = 0.0; y4[i] = 0.0;
        }
        csr_matvec(csr_p_size, stif3_csr_p, stif3_csr_j, stif3_csr_x,   W_now, y1); // dphidx * dphidy *   W_now
        csr_matvec(csr_p_size, stif4_csr_p, stif4_csr_j, stif4_csr_x,   W_now, y2); // dphidy * dphidx *   W_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, U2t_now, y3); // phi    * phi    * U1t_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x,  U2_now, y4); // phi    * phi    *  U2_now
        for (i = 0; i < node_num; i++)
        {
            rhs2[i]  =  - c[1][i] * y1[i] - c[3][i] *  y2[i] - rho[i] * mpml_dx[i] * y3[i] - rho[i] * mpml_dy[i] * y3[i] - rho[i] * mpml_dx[i] * mpml_dy[i] * y4[i]; 
        }

       /********************************************************************************************************************************************
        Equation 3:
         mass * U3tt_new = - c44 * dphidy * dphidy * U_now - 2.0 * rho * mpml_dy * phi * phi * U3t_now - rho * mpml_dy * mpml_dy * phi * phi * U3_now 
                           + phi * phi * Lx3_now + phi * phi * Lx4_now
        *********************************************************************************************************************************************/
        for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; y3[i] = 0.0; y4[i] = 0.0; y5[i] = 0.0;
        }
        csr_matvec(csr_p_size, stif2_csr_p, stif2_csr_j, stif2_csr_x,   U_now, y1); // dphidy * dphidy *  U_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, U3t_now, y2); // phi    * phi    * U1t_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x,  U3_now, y3); // phi    * phi    * U1_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Lx3_now, y4); // phi    * phi    * Lx1_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Lx4_now, y5); // phi    * phi    * Lx2_now
        for (i = 0; i < node_num; i++)
        {
            rhs3[i]  =  - c[3][i] * y1[i] - 2.0 * rho[i] * mpml_dy[i] * y2[i] -  rho[i] * mpml_dy[i] * mpml_dy[i] * y3[i] + y4[i] + y5[i];
        }

        /********************************************************************************************************************************************
         Equation 4:
         mass * Lx1_new = - dt * c11 * mpml_dxx * phi * dphidx * U_now - dt * rho * mpml_dx * phi * phi * Lx1_now - phi * phi * Lx1_now
        *********************************************************************************************************************************************/
        for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; 
        }
        csr_matvec(csr_p_size, stif5_csr_p, stif5_csr_j, stif5_csr_x,   U_now, y1); // phi    * dphidx *   U_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Lx1_now, y2); // phi    * phi    * Lx1_now
        for (i = 0; i < node_num; i++)
        {
            rhs4[i]  =  - dt * c[0][i] * mpml_dxx[i] * y1[i] - dt * rho[i] * mpml_dx[i] * y2[i] - y2[i];
        }

       /********************************************************************************************************************************************
         Equation 5:
         mass * Lx2_new = - dt * c44 * mpml_dyy_pxy * phi * dphidx * W_now - dt * rho * mpml_dy * phi * phi * Lx2_now - phi * phi * Lx2_now
        *********************************************************************************************************************************************/
        for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; 
        }
        csr_matvec(csr_p_size, stif5_csr_p, stif5_csr_j, stif5_csr_x,   W_now, y1); // phi    * dphidx *   W_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Lx2_now, y2); // phi    * phi    * Lx2_now
        for (i = 0; i < node_num; i++)
        {
            rhs5[i]  =  - dt * c[3][i] * mpml_dyy_pxy[i] * y1[i] - dt * rho[i] * mpml_dx[i] * y2[i] - y2[i];
        }
        
        /********************************************************************************************************************************************
         Equation 6:
         mass * Lx3_new = - dt * c13 * mpml_dxx_pyx * phi * dphidy * W_now - dt * rho * mpml_dx * phi * phi * Lx3_now - phi * phi * Lx3_now
        *********************************************************************************************************************************************/
        for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; 
        }
        csr_matvec(csr_p_size, stif6_csr_p, stif6_csr_j, stif6_csr_x,   W_now, y1); // phi    * dphidy *   W_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Lx3_now, y2); // phi    * phi    * Lx3_now
        for (i = 0; i < node_num; i++)
        {
            rhs6[i]  =  - dt * c[1][i] * mpml_dxx_pyx[i] * y1[i] - dt * rho[i] * mpml_dx[i] * y2[i] - y2[i];
        }

       /********************************************************************************************************************************************
         Equation 7:
         mass * Lx4_new = - dt * c44 * mpml_dyy * phi * dphidy * U_now - dt * rho * mpml_dy * phi * phi * Lx4_now - phi * phi * Lx4_now
        *********************************************************************************************************************************************/
        for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; 
        }
        csr_matvec(csr_p_size, stif6_csr_p, stif6_csr_j, stif6_csr_x,   U_now, y1); // phi    * dphidy *   U_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Lx4_now, y2); // phi    * phi    * Lx4_now
        for (i = 0; i < node_num; i++)
        {
            rhs7[i]  =  - dt * c[3][i] * mpml_dyy[i] * y1[i] - dt * rho[i] * mpml_dx[i] * y2[i] - y2[i];
        }