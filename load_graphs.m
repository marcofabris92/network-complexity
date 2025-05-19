% load graphs

clear all
close all
clc

load("graphs_ER_SW_BA_k4.mat","ER2","params_ER2");
ER2_k4 = ER2;
n = params_ER2(1);
ERthr_k4 = params_ER2(2);
                
load("graphs_ER_SW_BA_k4.mat","SW","params_SW");
SW_k4 = SW;
K_k4 = params_SW(2);
C_G_min_k4 = params_SW(3);
connSW_k4 = params_SW(4);

load("graphs_ER_SW_BA_k4.mat","BA","params_BA");
BA_k4 = BA;
n_ini_k4 = params_BA(1);
etaBA_k4 = params_BA(2);
m_ini_k4 = params_BA(3);
dnew_k4 = params_BA(5);
connBA_k4 = params_BA(6);





load("graphs_ER_SW_BA_k8.mat","ER2","params_ER2");
ER2_k8 = ER2;
ERthr_k8 = params_ER2(2);
                
load("graphs_ER_SW_BA_k8.mat","SW","params_SW");
SW_k8 = SW;
K_k8 = params_SW(2);
C_G_min_k8 = params_SW(3);
connSW_k8 = params_SW(4);

load("graphs_ER_SW_BA_k8.mat","BA","params_BA");
BA_k8 = BA;
n_ini_k8 = params_BA(1);
etaBA_k8 = params_BA(2);
m_ini_k8 = params_BA(3);
dnew_k8 = params_BA(5);
connBA_k8 = params_BA(6);


figure
hold on
subplot(2,3,1); plot(ER2_k4);
subplot(2,3,2); plot(SW_k4);
subplot(2,3,3); plot(BA_k4);
subplot(2,3,4); plot(ER2_k8);
subplot(2,3,5); plot(SW_k8);
subplot(2,3,6); plot(BA_k8);