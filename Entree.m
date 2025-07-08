function u = Entree(t,Mode)
global E
if Mode==1
    Si_in=25;
    Ss_in=80;
    Xi_in=40;
    Xs_in=200;
    Xbh_in=50;
    Xaob_in=0;
    Xnob_in=0;
    Xp_in=0;
    So_in=0;   
    SNo2_in=0;
    SNo3_in=0;
    SN2_in=0;
    SNh_in=30;
    SNd_in=5;
    XNd_in=10;
    Salk_in=4;
    u=[Si_in,Ss_in,Xi_in,Xs_in,Xbh_in,Xaob_in,Xnob_in,Xp_in,So_in,SNo2_in,SNo3_in,SN2_in,SNh_in,SNd_in,XNd_in,Salk_in];
elseif Mode==2
idx = find(abs(t-E(:,1))<10^(-2), 1);
u=E(idx,2:17);
end