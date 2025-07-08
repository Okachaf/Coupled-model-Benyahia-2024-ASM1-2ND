clc
close all
%%%%%%%%%%%%%%%%%%%%%%% MBR fouling model (Benyahia, 2024) (DOI: 10.3390/membranes14030069)
%%%%%%%%%%%%%%%%%%%%%%% coupled with ASM1-2ND model (PhD thesis by Farouk Aichouche) (HAL Id: tel-03986660)
%%%%%%%%%%%%%%%%%%%%%%% Okacha OUNADJELA
global  Volume Y_H Y_AOB Y_NOB f_P i_XB i_XP mu_H mu_AOB mu_NOB
global  b_H b_AOB b_NOB  K_S K_OH K_OA K_NO K_NH ny_g SOSAT k_a k_h K_X ny_NO2 ny_NO3
global  A0 sigma1 sigma2 alpha1 alpha2 ep TMP mu R0 Cs Cx beta
global  W1 W2 Csp Cxp McMax MpMax E
load dryinfluent.mat DRYINFLUENT
E=DRYINFLUENT;
N=60;%number of cycles
tf=[0 6.5/24];%filtration time
tb=[0 10/1440];%backwash time
Mode=1;% 1:const-inflow 2:dynamic-inflow 
%initialisation
SI_0 = 2;
SS_0 = 65;
XI_0 = 25;
XS_0 = 15;
XBH_0 = 10;
XAOB_0 = 7;
XNOB_0 = 3;
XP_0 = 0;
SO_0 = 0;
SNO2_0 = 0;
SNO3_0 = 0;
SN2_0 = 0;
SNH_0 = 0;
SND_0 = 0;
XND_0 = 1;
SALK_0 = 5;
MC_0=0;
MP_0=0;
C0=[SI_0,SS_0,XI_0,XS_0,XBH_0,XAOB_0,XNOB_0,XP_0,SO_0,SNO2_0,SNO3_0,SN2_0,SNH_0,SND_0,XND_0,SALK_0,MC_0,MP_0];
%inputs
Qw=.2;
KLA=250;
u=[Qw,KLA];
%fouling parameters
A0=32.4;%m2
sigma1=200;
sigma2=150;
alpha1=5*10^14*10^(-3);%(m/kg) to (m/g)
alpha2=10^13*10^(-3);%(m/kg) to (m/g)
ep=0.7;
TMP=.2;% bar
mu=0.001*10^(-5)/86400;%(Pa·s) to (bar·j)
R0=5.0897e+12;%(1/m)
Cs=0.02;
Cx=0.012;
beta=.015;
%backwash parameters
W1=25;
W2=5;
Csp=0.0001;
Cxp=0.001;
McMax=0;
MpMax=0;
%system parameters
Volume=3;%m3
mu_H = 4.0;  %6.0;
K_S = 10;  %20;
K_OH = 0.20;
K_NO = 0.5;
b_H = 0.3;  %0.62;
mu_AOB = 0.85;  %0.8;
mu_NOB = 0.65;
K_NH = 1;
K_OA = 0.4;
b_AOB = 0.17;
b_NOB = 0.15;%0.2;
ny_g = 0.8;
k_a = 0.05;  %0.08;
k_h = 3.0;
K_X = 0.1;  %0.03;
ny_NO2 = 0.8;  %0.4;
ny_NO3 = 0.38; 
Y_H = 0.45;
Y_AOB = 0.17;
Y_NOB = 0.15;
f_P = 0.08;
i_XB = 0.086;  %0.086;
i_XP = 0.06;
SOSAT=8;
%TSS
X_I2TSS = 0.75;
X_S2TSS = 0.75;
X_BH2TSS = 0.75;
X_AOB2TSS = 0.75;
X_NOB2TSS = 0.75;
X_P2TSS = 0.75;
%%%%%filtration 
[Te,Xe]=ode15s(@(tf,xe)filtration(tf,xe,u,Mode),tf,C0);
        T=Te;
        x=Xe;
        C0=Xe(end,:);
        Mc=Xe(:,17);
        Mp=Xe(:,18);
        A=A0./(1+(Mc/sigma1)+(Mp/sigma2));
        R=((alpha1*Mc./A)+(alpha2*Mp./(ep*A)));
        J=TMP./(mu*(R+R0));
        Qout=J.*A;
        

%%%%%filtration back-wash
  for n=2:2*N 
        
    if mod(n,2)==1
        
        t=tf+T(end)*ones(1,length(tf));
        [Te,Xe]=ode45(@(t,xe)filtration(t,xe,u,Mode),t,C0);
        T=[T;Te];
        C0=Xe(end,:);
        x=[x;Xe];

        Mc=Xe(:,17);
        Mp=Xe(:,18);
        
        AT=A0./(1+(Mc/sigma1)+(Mp/sigma2));
        RT=((alpha1*Mc./AT)+(alpha2*Mp./(ep*AT)));
        J=TMP./(mu*(RT+R0));
        A=[A;AT];
        R=[R;RT];
        Qout=[Qout;J.*AT];
        
        
    else 
        
        t=tb+T(end)*ones(1,length(tb));
        [Te,Xe]=ode45(@(t,xe)Back_wash(t,xe,u,Mode),t,C0);
        C0=Xe(end,:);
        T=[T;Te];
        x=[x;Xe];
        Mc=Xe(:,17);
        Mp=Xe(:,18);
        AT=A0./(1+(Mc/sigma1)+(Mp/sigma2));
        RT=((alpha1*Mc./AT)+(alpha2*Mp./(ep*AT)));
        J=0*TMP./(mu*(RT+R0));
        A=[A;AT];
        R=[R;RT];
        Qout=[Qout;J.*AT];
        
       
    end
  end
Qin=Qout+Qw*ones(length(Qout),1);

TSS = X_I2TSS * x(:,2) + X_S2TSS * x(:,3) + X_BH2TSS * x(:,4) + X_AOB2TSS * x(:,5) + X_NOB2TSS * x(:,6) + X_P2TSS * x(:,7);
Si = x(:,1);
Ss = x(:,2);
Xi = x(:,3);
Xs = x(:,4);
Xbh = x(:,5);
Xaob = x(:,6);
Xnob = x(:,7);
Xp = x(:,8);
So = x(:,9);   
SNo2 = x(:,10);
SNo3 = x(:,11);
SN2 = x(:,12);
SNh = x(:,13);
SNd = x(:,14);
XNd = x(:,15);
Salk = x(:,16);
Mc=x(:,17);
Mp=x(:,18);

St=Si + Ss + SNo2 + SNo3  + SNh + SNd +So+SN2+Salk;
Xt=Xi + Xs + Xbh + Xaob + Xnob + Xp + XNd;

figure
plot(T, St, 'b', 'LineWidth', 1.5)
grid on
title('St (g(COD)/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, Xt, 'b', 'LineWidth', 1.5)
grid on
title('xt (g(COD)/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, Si, 'b', 'LineWidth', 1.5)
grid on
title('Si (g(COD)/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, Ss, 'b', 'LineWidth', 1.5)
grid on
title('Ss (g(COD)/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, Xi, 'b', 'LineWidth', 1.5)
grid on
title('Xi (g(COD)/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, Xs, 'b', 'LineWidth', 1.5)
grid on
title('Xs (g(COD)/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, Xbh, 'b', 'LineWidth', 1.5)
grid on
title('Xbh (g(COD)/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, Xaob, 'b', 'LineWidth', 1.5)
grid on
title('Xaob (g(COD)/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, Xnob, 'b', 'LineWidth', 1.5)
grid on
title('Xnob (g(COD)/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, Xp, 'b', 'LineWidth', 1.5)
grid on
title('Xp (g(COD)/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, So, 'b', 'LineWidth', 1.5)
grid on
title('So (g(COD)/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, SNo2, 'b', 'LineWidth', 1.5)
grid on
title('SNo_2 (gN/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, SNo3, 'b', 'LineWidth', 1.5)
grid on
title('SNo_3 (gN/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, SN2, 'b', 'LineWidth', 1.5)
grid on
title('SN_2 (gN/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, SNh, 'b', 'LineWidth', 1.5)
grid on
title('SN_h (gN/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, SNd, 'b', 'LineWidth', 1.5)
grid on
title('SN_d (gN/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, XNd, 'b', 'LineWidth', 1.5)
grid on
title('XNd (gN/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, Salk, 'b', 'LineWidth', 1.5)
grid on
title('Salk (g(COD)/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, TSS, 'b', 'LineWidth', 1.5)
grid on
title('TSS (g(COD)/m^3)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, Qout, 'b', 'LineWidth', 1.5)
grid on
title('Qout (m^3/j)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)
ylabel('Qout (m^3/j)', 'FontSize', 12)

figure
plot(T, Qin, 'b', 'LineWidth', 1.5)
grid on
title('Qin (m^3/j)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, Mc, 'b', 'LineWidth', 1.5)
grid on
title('Mc (g)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, Mp, 'b', 'LineWidth', 1.5)
grid on
title('Mp (g)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, A, 'b', 'LineWidth', 1.5)
grid on
title('A (m^2)', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)

figure
plot(T, R, 'b', 'LineWidth', 1.5)
grid on
title('R', 'FontSize', 14)
xlabel('Temps (jours)', 'FontSize', 12)


%close all