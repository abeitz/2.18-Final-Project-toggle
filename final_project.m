%% Dynamic
clear all
param1 = [0.25 10 1 1 3 1 5 10 5 1 1 10 10];
[t1,u1pro1] = ode45(@(t,u) protein_ode(u,param1),[0 100],[0 0 0 200 200 0])
[t2,u1pro2] = ode45(@(t,u) protein_ode(u,param1),[0 100],[0 0 200 0 0 200])
figure(1)
clf
subplot(1,2,1)
plot(t1,u1pro1)
legend()
subplot(1,2,2)
plot(t2,u1pro2)
legend()

%% Steady State Analysis
clear all
param1 = [0.25 2 1 1 3 1 5 10 5 1 1 10 10];
u01m = [0 0 0 0 0 0 100 0];
u02m = [0 0 0 0 0 0 0 100];
u01p = [0 0 0 200 200 0];
u02p = [0 0 200 0 0 200];

param = [0.25 2 1 1 3 1 5 10 5 1 1 10 10];
al_vec = [1 2 3 4 5 6 7 8 9 10];
fig = figure(2)
clf
subplot(1,2,2)
set(gca, 'YScale', 'log','XScale','log')
hold on
for i = 1: length(al_vec)
    param(2) = al_vec(i);
    [t,u1miR] = ode45(@(t,u) miR_ode(u,param),[0 100],u01m);
    miR(2*i-1,1) = u1miR(end,7);
    miR(2*i-1,2) = u1miR(end,8);
    miR(2*i-1,3) = al_vec(i);
    [t,u2miR] = ode45(@(t,u) miR_ode(u,param),[0 100],u02m);
    miR(2*i,1) = u2miR(end,7);
    miR(2*i,2) = u2miR(end,8);
    miR(2*i,3) = al_vec(i);
    
end
scatter(miR(:,1),miR(:,2),40,miR(:,3),'filled')
a = colorbar();
a.Label.String = '\alpha';
a.Label.FontSize = 15;
title('miRNA Repression')
xlabel('A1')
ylabel('A2')
xlim([1e-3 1e3])
ylim([1e-3 1e3])
%legend('Location','southwest')
subplot(1,2,1)
set(gca, 'YScale', 'log','XScale','log')
hold on
for i = 1: length(al_vec)
    param(2) = al_vec(i);
    [t,u1pro] = ode45(@(t,u) protein_ode(u,param),[0 100],u01p);
    pro(2*i-1,1) = u1pro(end,5);
    pro(2*i-1,2) = u1pro(end,6);
    pro(2*i-1,3) = al_vec(i);
    [t,u2pro] = ode45(@(t,u) protein_ode(u,param),[0 100],u02p);
    pro(2*i,1) = u2pro(end,5);
    pro(2*i,2) = u2pro(end,6);
    pro(2*i,3) = al_vec(i);
end
scatter(pro(:,1),pro(:,2),40,miR(:,3),'filled')
colorbar();
title('Transcriptional Repression')
xlabel('A1')
ylabel('A2')
xlim([1e-3 1e3])
ylim([1e-3 1e3])

%%
clear miR pro
param = [0.25 2 1 1 3 1 5 10 5 1 1 10 10];
del_vec = [0.1 0.5 1 2 3 4 5 6 7 8 9];
figure(3)
clf
subplot(1,2,2)
set(gca, 'YScale', 'log','XScale','log')
hold on
for i = 1: length(del_vec)
    param(10) = del_vec(i);
    [t,u1miR] = ode45(@(t,u) miR_ode(u,param),[0 100],u01m);
    miR(2*i-1,1) = u1miR(end,7);
    miR(2*i-1,2) = u1miR(end,8);
    miR(2*i-1,3) = del_vec(i);
    [t,u2miR] = ode45(@(t,u) miR_ode(u,param),[0 100],u02m);
    miR(2*i,1) = u2miR(end,7);
    miR(2*i,2) = u2miR(end,8);
    miR(2*i,3) = del_vec(i);
end
scatter(miR(:,1),miR(:,2),40,miR(:,3),'filled')
a = colorbar();
a.Label.String = '\delta';
a.Label.FontSize = 15;
title('miRNA Repression')
xlabel('A1')
ylabel('A2')
xlim([1e-3 1e3])
ylim([1e-3 1e3])
subplot(1,2,1)
set(gca, 'YScale', 'log','XScale','log')
hold on
for i = 1: length(del_vec)
    param(10) = del_vec(i);
    [t,u1pro] = ode45(@(t,u) protein_ode(u,param),[0 100],u01p);
    pro(2*i-1,1) = u1pro(end,5);
    pro(2*i-1,2) = u1pro(end,6);
    pro(2*i-1,3) = del_vec(i);
    [t,u2pro] = ode45(@(t,u) protein_ode(u,param),[0 100],u02p);
    pro(2*i,1) = u2pro(end,5);
    pro(2*i,2) = u2pro(end,6);
    pro(2*i,3) = del_vec(i);
end
scatter(pro(:,1),pro(:,2),40,miR(:,3),'filled')
colorbar();
title('Transcriptional Repression')
xlabel('A1')
ylabel('A2')
xlim([1e-3 1e3])
ylim([1e-3 1e3])

%%
clear miR pro
param = [0.25 2 1 1 3 1 5 10 5 1 1 10 10];
gam_vec = [0.1 0.5 1 2 3 4 5 6 7 8 9];
figure(4)
clf
subplot(1,2,2)
set(gca, 'YScale', 'log','XScale','log')
hold on
for i = 1: length(gam_vec)
    param(11) = gam_vec(i);
    [t,u1miR] = ode45(@(t,u) miR_ode(u,param),[0 100],u01m);
    miR(2*i-1,1) = u1miR(end,7);
    miR(2*i-1,2) = u1miR(end,8);
    miR(2*i-1,3) = gam_vec(i);
    [t,u2miR] = ode45(@(t,u) miR_ode(u,param),[0 100],u02m);
    miR(2*i,1) = u2miR(end,7);
    miR(2*i,2) = u2miR(end,8);
    miR(2*i,3) = gam_vec(i);
end
scatter(miR(:,1),miR(:,2),40,miR(:,3),'filled')
a = colorbar();
a.Label.String = '\gamma';
a.Label.FontSize = 15;
title('miRNA Repression')
xlabel('A1')
ylabel('A2')
xlim([1e-3 1e3])
ylim([1e-3 1e3])
subplot(1,2,1)
set(gca, 'YScale', 'log','XScale','log')
hold on
for i = 1: length(gam_vec)
    param(11) = gam_vec(i);
    [t,u1pro] = ode45(@(t,u) protein_ode(u,param),[0 100],u01p);
    pro(2*i-1,1) = u1pro(end,5);
    pro(2*i-1,2) = u1pro(end,6);
    pro(2*i-1,3) = gam_vec(i);
    [t,u2pro] = ode45(@(t,u) protein_ode(u,param),[0 100],u02p);
    pro(2*i,1) = u2pro(end,5);
    pro(2*i,2) = u2pro(end,6);
    pro(2*i,3) = gam_vec(i);
end
scatter(pro(:,1),pro(:,2),40,miR(:,3),'filled')
colorbar();
title('Transcriptional Repression')
xlabel('A1')
ylabel('A2')
xlim([1e-3 1e3])
ylim([1e-3 1e3])

%%
clear miR pro
param = [0.25 2 1 1 3 1 5 10 5 1 1 10 10];
copy_vec = [1 5 10 50 100];
figure(5)
clf
subplot(1,2,2)
set(gca, 'YScale', 'log','XScale','log')
hold on
for i = 1: length(copy_vec)
    param(12) = copy_vec(i);
    param(13) = copy_vec(i);
    [t,u1miR] = ode45(@(t,u) miR_ode(u,param),[0 100],u01m);
    miR(2*i-1,1) = u1miR(end,7);
    miR(2*i-1,2) = u1miR(end,8);
    miR(2*i-1,3) = copy_vec(i);
    [t,u2miR] = ode45(@(t,u) miR_ode(u,param),[0 100],u02m);
    miR(2*i,1) = u2miR(end,7);
    miR(2*i,2) = u2miR(end,8);
    miR(2*i,3) = copy_vec(i);
    miR
end
scatter(miR(:,1),miR(:,2),40,miR(:,3),'filled')
a = colorbar();
a.Label.String = 'Copy Number';
a.Label.FontSize = 15;
title('miRNA Repression')
xlabel('A1')
ylabel('A2')
xlim([1e-3 1e3])
ylim([1e-3 1e3])
subplot(1,2,1)
set(gca, 'YScale', 'log','XScale','log')
hold on
for i = 1: length(copy_vec)
    param(12) = copy_vec(i);
    param(13) = copy_vec(i);
    [t,u1pro] = ode45(@(t,u) protein_ode(u,param),[0 100],u01p);
    pro(2*i-1,1) = u1pro(end,5);
    pro(2*i-1,2) = u1pro(end,6);
    pro(2*i-1,3) = copy_vec(i);
    [t,u2pro] = ode45(@(t,u) protein_ode(u,param),[0 100],u02p);
    pro(2*i,1) = u2pro(end,5);
    pro(2*i,2) = u2pro(end,6);
    pro(2*i,3) = copy_vec(i);
    pro
end
scatter(pro(:,1),pro(:,2),40,miR(:,3),'filled')
colorbar();
title('Transcriptional Repression')
xlabel('A1')
ylabel('A2')
xlim([1e-3 1e3])
ylim([1e-3 1e3])

%%
clear miR pro
param = [0.25 2 1 1 3 1 5 10 5 1 1 10 10];
splice_vec = [1 3 5 7 9 11 13 15 17 19];
figure(6)
clf
set(gca, 'YScale', 'log','XScale','log')
hold on
for i = 1: length(splice_vec)
    param(7) = splice_vec(i);
    [t,u1miR] = ode45(@(t,u) miR_ode(u,param),[0 100],u01m);
    miR(2*i-1,1) = u1miR(end,7);
    miR(2*i-1,2) = u1miR(end,8);
    miR(2*i-1,3) = splice_vec(i);
    [t,u2miR] = ode45(@(t,u) miR_ode(u,param),[0 100],u02m);
    miR(2*i,1) = u2miR(end,7);
    miR(2*i,2) = u2miR(end,8);
    miR(2*i,3) = splice_vec(i);
    miR
end
scatter(miR(:,1),miR(:,2),40,miR(:,3),'filled')
a = colorbar();
a.Label.String = 'k_s';
a.Label.FontSize = 15;
title('miRNA Repression')
xlabel('A1')
ylabel('A2')
xlim([1e-3 1e3])
ylim([1e-3 1e3])


%%
clear miR pro
param = [0.25 2 1 1 3 1 5 10 5 1 1 10 10];
ap_vec = [1 3 5 7 9 11 13 15 17 19];
figure(7)
clf
set(gca, 'YScale', 'log','XScale','log')
hold on
for i = 1: length(ap_vec)
    param(5) = ap_vec(i);
    [t,u1miR] = ode45(@(t,u) miR_ode(u,param),[0 100],u01m);
    miR(2*i-1,1) = u1miR(end,7);
    miR(2*i-1,2) = u1miR(end,8);
    miR(2*i-1,3) = ap_vec(i);
    [t,u2miR] = ode45(@(t,u) miR_ode(u,param),[0 100],u02m);
    miR(2*i,1) = u2miR(end,7);
    miR(2*i,2) = u2miR(end,8);
    miR(2*i,3) = ap_vec(i);
    miR
end
scatter(miR(:,1),miR(:,2),40,miR(:,3),'filled')
a = colorbar();
a.Label.String = 'miRNA K_d';
a.Label.FontSize = 15;
title('miRNA Repression')
xlabel('A1')
ylabel('A2')
xlim([1e-3 1e3])
ylim([1e-3 1e3])

%% Frequency Response
clear all
p = [0.25 4 10 10 30 10 4 10 5 10 10 10 10];
u01m = [0 0 0 0 0 0 100 0];
u01p = [0 0 0 200 200 0];
[t,u1miR] = ode45(@(t,u) miR_ode(u,p),[0 100],u01m);
[t,u1pro] = ode45(@(t,u) protein_ode(u,p),[0 100],u01p);
miRE = u1miR(end,:);
proE = u1pro(end,:);
Am = [-(p(7)+p(10))-(p(8)/(p(6)+p(8)))*p(5)*miRE(6), 0, 0, 0, 0, -(p(8)/(p(6)+p(8)))*p(5)*miRE(1), (p(2)-p(1))*(p(3)*p(4)*p(12)/(p(3)*miRE(7)+p(4))^2), 0;
    0, -(p(7)+p(10))-(p(8)/(p(6)+p(8)))*p(5)*miRE(5), 0, 0, -(p(8)/(p(6)+p(8)))*p(5)*miRE(2), 0, 0, (p(2)-p(1))*(p(3)*p(4)*p(13)/(p(3)*miRE(8)+p(4))^2);
    p(7), 0, -(p(8)/(p(6)+p(8)))*p(5)*miRE(6)-p(10), 0, 0, -(p(8)/(p(6)+p(8)))*p(5)*miRE(3), 0, 0;
    0, p(7), 0, -(p(8)/(p(6)+p(8)))*p(5)*miRE(5)-p(10), -(p(8)/(p(6)+p(8)))*p(5)*miRE(4), 0, 0, 0;
    p(7), 0, 0, 0, -p(10), 0, 0, 0;
    0, p(7), 0, 0, 0, -p(10), 0, 0;
    0, 0, p(9), 0, 0, 0, -p(11), 0;
    0, 0, 0, p(9), 0, 0, 0, -p(11)]
Bm_txn = [1; 1; 0; 0; 0; 0; 0; 0];
Bm_tln = [0; 0; 0; 0; 0; 0; 1; 1];
Cm = [0, 0, 0, 0, 0, 0, 1, 0];
Ap = [-p(10), 0, (p(1)+(p(2)*p(3)*proE(5)/p(4)))*(p(12)*p(3)/p(4))/(1+(p(3)*proE(5)/p(4))+(p(3)*proE(3)/p(4)))^2, 0, p(12)*((p(2)*p(3)/p(4))*(1+(p(3)*proE(3)/p(4)))-(p(1)*p(3)/p(4)))/(1+(p(3)*proE(5)/p(4))+(p(3)*proE(3)/p(4)))^2, 0;
    0, -p(10), 0, (p(1)+(p(2)*p(3)*proE(6)/p(4)))*(p(13)*p(3)/p(4))/(1+(p(3)*proE(6)/p(4))+(p(3)*proE(4)/p(4)))^2, 0,  p(13)*((p(2)*p(3)/p(4))*(1+(p(3)*proE(4)/p(4)))-(p(1)*p(3)/p(4)))/(1+(p(3)*proE(6)/p(4))+(p(3)*proE(4)/p(4)))^2;
    0, p(9), -p(11), 0, 0, 0;
    p(9), 0, 0, -p(11), 0, 0;
    p(9), 0, 0, 0, -p(11), 0;
    0, p(9), 0, 0, 0, -p(11)]
Bp_txn = [1; 1; 0; 0; 0; 0];
Bp_tln = [0; 0; 1; 1; 1; 1];
Cp = [0, 0, 0, 0, 1, 0];

%%
s = tf('s');
figure(8)
subplot(1,2,1)
Mm_txn = Cm*inv(s*eye(8)-Am)*Bm_txn
Mp_txn = Cp*inv(s*eye(6)-Ap)*Bp_txn
bode(Mp_txn,Mm_txn)
legend('Design 1','Design 2')
title('Transcriptional Noise')
subplot(1,2,2)
Mm_tln = Cm*inv(s*eye(8)-Am)*Bm_tln
Mp_tln = Cp*inv(s*eye(6)-Ap)*Bp_tln
bode(Mp_tln,Mm_tln)
legend('Design 1','Design 2')
title('Translational Noise')
%%
figure(10)
clf
margin(Mm_txn)
%% Helper Functions
function F = miR_ode(u,p)
% u = [t1 t2 m1 m2 r1 r2 A1 A2]

% Define parameters
al0 = p(1); al= p(2); a = p(3); d = p(4);
ap = p(5); dp = p(6); ks = p(7); km = p(8);
kap = p(9); del = p(10); gam = p(11); D1 = p(12); D2 = p(13);

% Define ODEs
F = [al0*D1 + (al-al0)*(a*u(7)*D1)/(a*u(7)+d) - (ks+del)*u(1) - (km/(dp+km))*ap*u(1)*u(6);
    al0*D2 + (al-al0)*(a*u(8)*D2)/(a*u(8)+d) - (ks+del)*u(2) - (km/(dp+km))*ap*u(2)*u(5);
    ks*u(1) - (km/(dp+km))*ap*u(3)*u(6) - del*u(3);
    ks*u(2) - (km/(dp+km))*ap*u(4)*u(5) - del*u(4);
    ks*u(1) - del*u(5);
    ks*u(2) - del*u(6);
    kap*u(3) - gam*u(7);
    kap*u(4) - gam*u(8)];
end

function F = protein_ode(u,p)
% u = [m1 m2 R1 R2 A1 A2]

% Define parameters
al0 = p(1); al= p(2); a = p(3); d = p(4);
ap = p(5); dp = p(6); ks = p(7); km = p(8);
kap = p(9); del = p(10); gam = p(11); D1 = p(12); D2 = p(13);

% Define ODEs
F = [(al0 + al*(a*u(5)/d))*(D1/(1+(a*u(5)/d)+(a*u(3)/d))) - del*u(1);
    (al0 + al*(a*u(6)/d))*(D1/(1+(a*u(6)/d)+(a*u(4)/d))) - del*u(2);
    kap*u(2) - gam*u(3);
    kap*u(1) - gam*u(4);
    kap*u(1) - gam*u(5);
    kap*u(2) - gam*u(6)];
end
