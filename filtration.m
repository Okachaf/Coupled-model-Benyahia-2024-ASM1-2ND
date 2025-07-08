function sys = filtration(t,x,u,Mode)
global  Volume Y_H Y_AOB Y_NOB f_P i_XB i_XP mu_H mu_AOB mu_NOB
global  b_H b_AOB b_NOB  K_S K_OH K_OA K_NO K_NH ny_g SOSAT k_a k_h K_X ny_NO2 ny_NO3
global  A0 sigma1 sigma2 alpha1 alpha2 ep TMP mu R0 Cs Cx beta McMax MpMax 

       Si = x(1);
       Ss = x(2);
       Xi = x(3);
       Xs = x(4);
       Xbh = x(5);
       Xaob = x(6);
       Xnob = x(7);
       Xp = x(8);
       So = x(9);   
       SNo2 = x(10);
       SNo3 = x(11);
       SN2 = x(12);
       SNh = x(13);
       SNd = x(14);
       XNd = x(15);
       Salk = x(16);
       Mc=x(17);
       Mp=x(18);

       St=Si + Ss + SNo2 + SNo3  + SNh + SNd +So+SN2+ Salk;
       Xt=Xi + Xs + Xbh + Xaob + Xnob + Xp + XNd;

       usys=Entree(t,Mode);

       Si_in=usys(1);
       Ss_in=usys(2);
       Xi_in=usys(3);
       Xs_in=usys(4);
       Xbh_in=usys(5);
       Xaob_in=usys(6);
       Xnob_in=usys(7);
       Xp_in=usys(8);
       So_in=usys(9);   
       SNo2_in=usys(10);
       SNo3_in=usys(11);
       SN2_in=usys(12);
       SNh_in=usys(13);
       SNd_in=usys(14);
       XNd_in=usys(15);
       Salk_in=usys(16);
       Qw=u(1);
       KLA=u(2);
       


              proc1 = mu_H * (((Ss / (K_S + Ss)) * (So / (K_OH + So))))  * Xbh; 
proc2 = mu_H * ((Ss / (K_S + Ss)) * (K_OH / (K_OH + So)) * (SNo2 / (K_NO + SNo2)) * (SNh / (K_NH + SNh))) * ny_NO2 * Xbh;
proc3 = mu_H * ((Ss / (K_S + Ss)) * (K_OH / (K_OH + So)) * (SNo3 / (K_NO + SNo3)) * (SNh / (K_NH + SNh))) * ny_NO3 * Xbh;
proc4 = b_H * Xbh;
proc5 = mu_AOB * (((SNh) / (K_NH + SNh)) * ((So) / (K_OA + So))) * Xaob;
proc6 = mu_NOB * (((SNo2) / (K_NO + SNo2)) * ((So) / (K_OA + So)) * ((K_NH) / (K_NH + SNh))) * Xnob;
proc7 = b_AOB * Xaob;
proc8 = b_NOB * Xnob;
proc9 = k_a * SNd * Xbh;
proc10 = k_h * (((Xs / Xbh) / (K_X + (Xs / Xbh))) * ((So / (K_OH + So)) + ny_g * (K_OH / (K_OH + So)) * ((SNo2 + SNo3) / (K_NO + SNo2 + SNo3)))) * Xbh;
proc11 = proc10 * (XNd / Xs);
       react1 = 0.0;
  react2 = (-proc1 - proc2 - proc3) / Y_H + proc10;
  react3 = 0.0;
  react4 = (1.0 - f_P) * (proc4 + proc8 + proc7) - proc10;
  react5 =  ( proc1 + proc2 + proc3 )- proc4;
  react6 =proc5 - proc7;
  react7 = proc6 - proc8;
  react8 = f_P * (proc4 + proc7 + proc8);
  react9 = (-(1 - Y_H) / Y_H) * proc1 + ((-3.4268028 + Y_AOB) / Y_AOB) * proc5 + ((-1.1414285 +Y_NOB) / Y_NOB)*proc6;
  react10 = (-(1 - Y_H) / (1.7242857 * Y_H)) * proc2 + proc5 / Y_AOB -proc6 / Y_NOB;
  react11 = -((1 - Y_H) / (2.8571 * Y_H)) * proc3 + proc6 / Y_NOB;
  react12 = (1 - Y_H) / (1.7242857 * Y_H) * proc2 + ((1 - Y_H) / (2.8571 * Y_H)) * proc3;
  react13 = (-i_XB)* (proc1 + proc2 + proc3)-i_XB*proc6 + (-i_XB - (1.0 / Y_AOB)) * proc5 + proc9;
  react14 = -proc9 + proc11;
  react15 = (i_XB - f_P * i_XP) * (proc4 + proc7 + proc8) - proc11;
  react16 = (-i_XB * (1 / 14)) * (proc1) + (proc2 + proc3) * ((1 - Y_H) / (40* Y_H) - ( i_XB / 14))...
        + proc5 * (-(1 / 7 * Y_AOB) - (i_XB / 14)) - (i_XB / 14)*proc6 +(proc9 / 14);
           
           A=A0/(1+(Mc/sigma1)+(Mp/sigma2));
           R=((alpha1*Mc/A)+(alpha2*Mp/(ep*A)));
           J=TMP/(mu*(R+R0));
           Qout=J*A;
           Qin=Qout+Qw;


           dx(1)=(1/Volume)*(Si_in*Qin-Si*Qout-Si*Qw-Si*Cs*Qout)+react1;
           dx(2)=(1/Volume)*(Ss_in*Qin-Ss*Qout-Ss*Qw-Ss*Cs*Qout)+react2;
           dx(3)=(1/Volume)*(Xi_in*Qin-Xi*Qw-Xi*Cx*Qout)+react3;
           dx(4)=(1/Volume)*(Xs_in*Qin-Xs*Qw-Xs*Cx*Qout)+react4;
           dx(5)=(1/Volume)*(Xbh_in*Qin-Xbh*Qw-Xbh*Cx*Qout)+react5;
           dx(6)=(1/Volume)*(Xaob_in*Qin-Xaob*Qw-Xaob*Cx*Qout)+react6;
           dx(7)=(1/Volume)*(Xnob_in*Qin-Xnob*Qw-Xnob*Cx*Qout)+react7;
           dx(8)=(1/Volume)*(Xp_in*Qin-Xp*Qw-Xp*Cx*Qout)+react8;
           dx(9)=(1/Volume)*(So_in*Qin-So*Qout-So*Qw-So*Cs*Qout)+react9+KLA*(SOSAT-So);
           dx(10)=(1/Volume)*(SNo2_in*Qin-SNo2*Qout-SNo2*Qw-SNo2*Cs*Qout)+react10;
           dx(11)=(1/Volume)*(SNo3_in*Qin-SNo3*Qout-SNo3*Qw-SNo3*Cs*Qout)+react11;
           dx(12)=(1/Volume)*(SN2_in*Qin-SN2*Qout-SN2*Qw-SN2*Cs*Qout)+react12;
           dx(13)=(1/Volume)*(SNh_in*Qin-SNh*Qout-SNh*Qw-SNh*Cs*Qout)+react13;
           dx(14)=(1/Volume)*(SNd_in*Qin-SNd*Qout-SNd*Qw-SNd*Cs*Qout)+react14;
           dx(15)=(1/Volume)*(XNd_in*Qin-XNd*Qw-XNd*Cx*Qout)+react15;
           dx(16)=(1/Volume)*(Salk_in*Qin-Salk*Qout-Salk*Qw-Salk*Cs*Qout)+react16;
           dx(17)=Qout*(Cs*St+Cx*Xt);
           dx(18)=Qout*(beta*St);

        if Mc>McMax
            McMax=Mc;
        end
        if Mp>MpMax
            MpMax=Mp;
        end
        
sys=dx';
end