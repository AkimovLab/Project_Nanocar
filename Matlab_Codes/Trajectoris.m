% This code plots the trajectories relating to motion of C60 at different temperatures
clc, clear
load allT4.mat

% Data column order:
%1Step 2KinEng 3PotEng 4Temp 5Press 6c_ke_nc 7c_lennard 8c_pe_nc 9v_xcm_x 10v_xcm_y 11v_xcm_z 
%12v_vcm_x v_vcm_y v_vcm_z v_wcm_x v_wcm_y v_wcm_z v_fcm_x v_fcm_y v_fcm_z 
%21c_msdd[1] c_msdd[2] c_msdd[3] v_x1_x v_x1_y v_x1_z v_x2_x v_x2_y v_x2_z 
%30v_wheel1_x v_wheel1_y v_wheel1_z v_wheel1_vx v_wheel1_vy v_wheel1_vz v_wheel1_ox v_wheel1_oy v_wheel1_oz 
%39v_wheel2_x v_wheel2_y v_wheel2_z v_wheel2_vx v_wheel2_vy v_wheel2_vz v_wheel2_ox v_wheel2_oy v_wheel2_oz 
%48v_wheel3_x v_wheel3_y v_wheel3_z v_wheel3_vx v_wheel3_vy v_wheel3_vz v_wheel3_ox v_wheel3_oy v_wheel3_oz 
%57v_wheel4_x v_wheel4_y v_wheel4_z v_wheel4_vx v_wheel4_vy v_wheel4_vz v_wheel4_ox v_wheel4_oy v_wheel4_oz 
%66v_chassi_x v_chassi_y v_chassi_z v_chassi_vx v_chassi_vy v_chassi_vz v_chassi_ox v_chassi_oy v_chassi_oz 
%75c_temp_nc c_temp_sub 

% Data column order:
%1Step 2CPU 3PotEng 4KinEng 5Temp 6Lx 7Ly 8Press 
%9v_xc_x 10v_xc_y 11v_xc_z 12c_pe_c60 13c_lennard 14c_ke_c60 
%15v_vc_x 16v_vc_y 17v_vc_z
%18v_x1_x 19v_x1_y 20v_x1_z 21v_x2_x 22v_x2_y 23v_x2_z 24c_pe_sub 25c_ke_sub 
%26v_wc_x 27v_wc_y 28v_wc_z 29v_w12_x 30v_w12_y 31v_w12_z 32c_temp_c60 33c_temp_sub


T=[5 10 30 50 75 100 150 200 300 400 500 600 700 800 900 1000];

s = 20000; % start point 10000

x(:,:)   =imdata(s:end,9,:);
y(:,:)   =imdata(s:end,10,:);
z(:,:)   =imdata(s:end,11,:);
len(:,:) =imdata(s:end,7,:);
pe(:,:)  =imdata(s:end,8,:);
potE(:,:)=imdata(s:end,3,:);
NP       = length(x);           % number of data points
thermo = 200;
t        = imdata(1:NP,1);
dt = 1e-3;                        % 0.001 ps time between trajectory points
time=dt*t;                        % time in (ps)

%% Trajectory
nsub=length(T);
col=jet(nsub); % parula, hsv, hot, pink, flag
               %  autumn, bone, colorcube, cool, copper, gray, 
               % jet, lines, prism, spring, summer, winter    
iall=1:16;
iva=16:-1:1;
i1=1:4;
i2=5:10;
i3=11:16;

xt=x;yt=y;

yt(:,3)=y(:,3)-10;
xt(:,4)=x(:,4)-25;
xt(:,1:2)=x(:,1:2)+30;
yt(:,4)=y(:,4)+30;
 
figure(1)
for p=i1
% Ekave(j)=mean(imdata(10:end,6,j));
% Etot=imdata(:,3,j)+imdata(:,6,j);
% Etot2(j)=mean(imdata(10:end,9,j));
% plot(imdata(:,1,j),Etot);
txt = ['T=',num2str(T(p)), 'K'];
plot(xt(:,p),yt(:,p),'LineWidth',3,'MarkerSize',12,'color',col(p,:),'DisplayName',txt); 
hold on
end
xlabel('x (\AA)','Interpreter','latex')
ylabel('y (\AA)','Interpreter','latex')
set(gca,'FontName','Cambria','FontSize',20);
legend ('Location','Best')

xt=x;yt=y;

xt(:,9)=x(:,9)-700;
yt(:,9)=y(:,9)-500;
xt(:,5:7)=x(:,5:7)-600;
xt(:,10)=x(:,10)-300;
yt(:,6)=y(:,6)+300;
yt(:,7)=y(:,7)+500;
xt(:,8)=x(:,8)-300;
yt(:,8)=y(:,8)+700;

figure(2)
for p=i2
txt = ['T=',num2str(T(p)), 'K'];
plot(xt(:,p),yt(:,p),'LineWidth',3,'MarkerSize',12,'color',col(p,:),'DisplayName',txt); 
hold on
end
xlabel('x (\AA)','Interpreter','latex')
ylabel('y (\AA)','Interpreter','latex')
legend ('Location','Best')
set(gca,'FontName','Cambria','FontSize',20);


xt=x;yt=y;

yt(:,13:14)=y(:,13:14)-2300;
xt(:,15)=x(:,15)-1500;
xt(:,14)=x(:,14)-2700;
yt(:,16)=y(:,16)+2000;
xt(:,16)=x(:,16)+600;
yt(:,15)=y(:,15)+500;
xt(:,13)=x(:,13)-500;
% xt(:,11)=x(:,11)+1200;
% yt(:,11)=y(:,11)-2000;
% xt(:,12)=x(:,12)-2000;
yt(:,12)=y(:,12)-1500;

figure(3)
for p=i3
txt = ['T=',num2str(T(p)), 'K'];
plot(xt(:,p),yt(:,p),'LineWidth',3,'MarkerSize',12,'color',col(p,:),'DisplayName',txt); 
hold on
end
xlabel('x (\AA)','Interpreter','latex')
ylabel('y (\AA)','Interpreter','latex')
set(gca,'FontName','Cambria','FontSize',20);
legend ('Location','Best')

%% x-t levy flight
figure(14)
for Ti=[5 9 12]
txt = ['T=',num2str(T(Ti)), 'K'];
plot(time/1000,x(:,Ti),'LineWidth',2,'MarkerSize',12,'color',col(Ti,:),'DisplayName',txt); 
%plot(time/1000,x(:,Ti),'.')
hold on
end
xlabel('t (ns)','Interpreter','latex')
ylabel('x (\AA)','Interpreter','latex')
set(gca,'FontName','Cambria','FontSize',20);
legend ('Location','Best')

% z-t
Ti=3;
figure(15)
plot(time-20000,z(:,Ti),'o','LineWidth',1,'MarkerSize',8); 
hold on
xlabel('t (ps)','Interpreter','latex')
ylabel('z (\AA)','Interpreter','latex')
% legend('T=30 K','T=300 K')
% legend('C60/SLG','C60/DLG','C60/FLG','Location','Best')
set(gca,'FontName','Cambria','FontSize',20);
axis([0,500,4.5,8])

Ti=9;
figure(16)
plot(time-20000,z(:,Ti),'o','LineWidth',1,'MarkerSize',8); 
hold on
xlabel('t (ps)','Interpreter','latex')
ylabel('z (\AA)','Interpreter','latex')
% legend('T=30 K','T=300 K')
% legend('NT/SLG','NT/DLG','NT/FLG','Location','Best')
%legend('NC/SLG','NC/DLG','NC/FLG','Location','Best')
set(gca,'FontName','Cambria','FontSize',20);
axis([0,500,4.5,8])

