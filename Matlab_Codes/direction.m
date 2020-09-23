clc, clear
load allT3.mat
T=[5 10 30 50 75 100 150 200 300 400 500 600 700 800 900 1000];
%  1 2  3  4  5  6   7   8   9   10  11  12  13  14  15  16

a = 20000; % start point 10000

vx(:,:)= imdata(a:end,12,:);      % throw off 5% of initial data points
vy(:,:)= imdata(a:end,13,:);      % vy(t,T)
wz(:,:)= imdata(a:end,17,:);
Nt = length(vx);                  % number of data points
t = imdata(1:Nt,1);
thermo=200;
dt = 1e-3;                        % 0.001 ps time between trajectory points
time=dt*t;                        % time in (ps)

% tgtetax = vy./vx;
% tetaatan=atan(tgtetax);
% vt=sqrt(vx.^2+vy.^2);
% costetax= vx./vt;
% tetaacos= acos(costetax);         % velocity and x 

x1(:,:)= imdata(a:end,24,:); 
y1(:,:)= imdata(a:end,25,:); 
x2(:,:)= imdata(a:end,27,:); 
y2(:,:)= imdata(a:end,28,:); 
rt= sqrt((x2-x1).^2+(y2-y1).^2);
% cosphi= (x2-x1)./rt;

x21= x2-x1;
y21= y2-y1;
% cosdev= (vx.*x21+vy.*y21)./(vt.*rt);

teta= atan2(vy, vx);              % motion direction
phi= atan2(y21, x21);             % chasis direction

% phiacos= acos(cosphi);               % chasis and x angle
% devacos= acos(cosdev);               % velocity and chasis
dev=phi-teta;

% figure(1)
% plot(teta(:,1),'o')
% hold on
% plot(tetaatan(:,1),'o')
% plot(tetaacos(:,1),'o')
% figure(2)
% plot(phi(:,5),'o')
% hold on
% plot(phiacos(:,5),'o')
% figure(3)
% plot(dev(:,5),'o')
% hold on
% plot(devacos(:,5),'o')

alfa = -pi:0.0175:pi;
alfad=-179:180;

i=[6 9 12]; %[1 5 9 12 16]; % T =  [1 6 9 12]
iall=1:16;
iva=[16 12 9 5 1];
nsub=length(T);
col=jet(nsub); % parula, hsv, hot, pink, flag
               %  autumn, bone, colorcube, cool, copper, gray, 
               % jet, lines, prism, spring, summer, white, winter, rgbplot, colormap 
               



%% population
for p=i    %i=[6 9 12] 100,300,600
    txt = ['T = ',num2str(T(p)), 'K'];
   
    figure(50)
    hphi= histogram(phi(:,p), 360);
    figure(51)
    hteta= histogram(teta(:,p), 360); 
    
    figure(100)
    polarplot(alfa,hteta.Values,'o-','LineWidth',2,'MarkerSize',2,'color',col(p,:),'DisplayName',txt)
    hold on
    
%     figure(p)
%     plot(hteta.Values,hphi.Values,'o','LineWidth',2,'MarkerSize',6,'color',col(p,:))
%     xlabel(' \theta population ','Interpreter','tex' )
%     ylabel(' \phi population ','Interpreter','tex' )
%     set(gca,'FontName','Cambria','FontSize',14)
%     tname = sprintf('T=%d K', T(p));
%     title(tname,'FontSize',24)
        
    figure(101)
    polarplot(alfa,hphi.Values,'o-','LineWidth',2,'MarkerSize',2,'color',col(p,:),'DisplayName',txt)
    hold on
    % ang_vec = angle (pi()*tetax(:,p)/180);
    % rose(pi()*tetax(:,p)/180 , 360)
end
figure(101)
title('chasis direction (\phi)','Interpreter','tex' )
set(gca,'FontName','Cambria','FontSize',16);
legend show

figure(100)
title('motion direction (\theta)','Interpreter','tex' )
set(gca,'FontName','Cambria','FontSize',16);
legend show


for p=11    %i=[6 9 12] 100,300,600
       
    figure(50)
    hphi= histogram(phi(:,p), 360);
    figure(51)
    hteta= histogram(teta(:,p), 360); 
    
    figure(75)
    polarplot(alfa,hteta.Values,'o-','LineWidth',2,'MarkerSize',2)
    hold on
        figure(76)
    polarplot(alfa,hphi.Values,'o-','LineWidth',2,'MarkerSize',2)
    hold on
    % ang_vec = angle (pi()*tetax(:,p)/180);
    % rose(pi()*tetax(:,p)/180 , 360)
end
figure(76)
title('chasis direction (\phi)','Interpreter','tex' )
set(gca,'FontName','Cambria','FontSize',16);
legend show

figure(75)
txt = ['motion direction (\theta)',' at T=',num2str(T(p)),' K'];
title(txt)
set(gca,'FontName','Cambria','FontSize',16);
legend show

% legend('NC/SLG','NC/DLG','NC/FLG')
% legend('NT/SLG','NT/DLG','NT/FLG')


%% deviation from chasis direction 
devm=dev;
% for k=iall
%     for j=1:Nt
%         if dev(j,k)>pi
%             devm(j,k)=2*pi-dev(j,k);
%         elseif dev(j,k)<-pi
%             devm(j,k)=2*pi+dev(j,k);
%         end
%     end
% end
for k=iall
    for j=1:Nt
        if dev(j,k)<0
            devm(j,k)=2*pi+dev(j,k);
        end
    end
end

    
for q=i
    figure(90)
    hdev= histogram(devm(:,q)*180/pi, 360);
    hold on
        
    txt = ['T = ',num2str(T(q)), 'K'];
    
    figure(91)
    polarplot(hdev.Values,'o-','LineWidth',2,'MarkerSize',2,'color',col(q,:),'DisplayName',txt)
    hold on
%     tname = sprintf('T=%d K', T(q));
%     legend(tname,'Location','Best')
    figure(92)
    rose(devm(:,q),360)
    hold on
end
% legend('T=5 K','T=75 K','T=300 K','T=600 K','T=1000 K','Location','Best')
title('deviance of motion from chasis direction (\theta-\phi)','Interpreter','tex' )
set(gca,'FontName','Cambria','FontSize',12);
legend show


%% wz

for p=iva 
    figure(80)
    plot(wz(:,p),teta(:,p),'.','LineWidth',2,'MarkerSize',9,'color',col(p,:))
    hold on
end
set(gca,'FontName','Cambria','FontSize',16);
legend('T=1000 K','T=600 K','T=300 K','T=75 K','T=55 K','Location','Best')
xlabel(' \omega z ','Interpreter','tex' )
ylabel('\theta (motion direction)','Interpreter','tex')

for p=iva 
    figure(81)
    plot(time/1000, wz(:,p),'.','LineWidth',2,'MarkerSize',9,'color',col(p,:))
    hold on
end
set(gca,'FontName','Cambria','FontSize',16);
legend('T=1000 K','T=600 K','T=300 K','T=75 K','T=55 K','Location','Best')
xlabel('time(ns)','Interpreter','tex')
ylabel(' \omega z ','Interpreter','tex' )


for q=iva
    txt = ['T = ',num2str(T(q)), 'K'];
    figure(72)
    plot(teta(:,q),phi(:,q),'.','LineWidth',2,'MarkerSize',9,'color',col(q,:),'DisplayName',txt)
    hold on

end
xlabel('motion direction (\theta)','Interpreter','tex')
ylabel('chasis direction (\phi)','Interpreter','tex' )
set(gca,'FontName','Cambria','FontSize',16);
legend show

%%
% costetax=zeros(Nt,length(T));
% costetay=zeros(Nt,length(T));
% costeta60=zeros(Nt,length(T));
% tgtetax=zeros(Nt,length(T));

% for i = 1:Nt;
%     for k=1:length(T)
%         costetax(i,k)=vx(i,k)./(sqrt(vx(i,k)^2+vy(i,k)^2));
%         costetay(i,k)=vy(i,k)./(sqrt(vx(i,k)^2+vy(i,k)^2));
%         costeta60(i,k)=(vx(i,k)+vy(i,k))./(sqrt(2*vx(i,k)^2+2*vy(i,k)^2));
%     end
% end

% tetay= acosd(costetay);
% teta60=acosd(costeta60);
% tetayx=90-tetax;
% sinx=sind(tetax);

% for i=[15 13 10 7 4 1] %length(T):-3:1
%     figure(90)
%     plot(time/1000,tetarv(:,i),'.','LineWidth',2,'MarkerSize',9)
%     hold on
% end
% set(gca,'FontName','Cambria','FontSize',16);
% legend('T=1000 K','T=700 K','T=400 K','T=150 K','T=50 K','T=5 K','Location','Best')
% xlabel('time(ns)','Interpreter','tex' )
% ylabel('( \theta ) velocity and chasis angle','Interpreter','tex')



% 
% for i=[15 13 10 7 4 1] %length(T):-3:1
%     figure(91)
%     plot(cosrx(:,i),costetax(:,i),'.','LineWidth',2,'MarkerSize',9) % wz
%     hold on
% end
% set(gca,'FontName','Cambria','FontSize',16);
% legend('T=1000 K','T=700 K','T=400 K','T=150 K','T=50 K','T=5 K','Location','Best')
% xlabel('cos( \theta2) chasis and x angle','Interpreter','tex' )
% ylabel('cos( \theta ) velocity and x angle','Interpreter','tex')

% for i=[15 13 10 7 4 1] %length(T):-3:1
%     figure(101)
%     plot(wz(:,i),tgtetax(:,i),'.','LineWidth',2,'MarkerSize',9)
%     hold on
% end
% set(gca,'FontName','Cambria','FontSize',16);
% legend('T=1000 K','T=700 K','T=400 K','T=150 K','T=50 K','T=5 K','Location','Best')
% xlabel(' \omega z ','Interpreter','tex' )
% ylabel('tg( \theta )','Interpreter','tex')


% for i=[15 13 10 7 4 1] %length(T):-3:1
%     figure(92)
%     plot(tetarx(:,i),tetax(:,i),'.','LineWidth',2,'MarkerSize',9)
%     hold on
% end
% set(gca,'FontName','Cambria','FontSize',16);
% legend('T=1000 K','T=700 K','T=400 K','T=150 K','T=50 K','T=5 K','Location','Best')
% xlabel('( \theta2) chasis and x angle','Interpreter','tex' )
% ylabel('( \theta ) velocity and x angle','Interpreter','tex')



% for p=1:16
%     figure(103)
%     hold on
%     plot(time(1:10000),costetax(1:10000,p),'.','LineWidth',2,'MarkerSize',9)
%     
%     figure(104)
%     hold on
%     histogram(costetax(:,p), 24)
% %     histogram(costetay(:,p), 24)
% %     histogram(costeta60(:,p), 24)
%     %plot(time(1:100:end),tetax(1:100:end,p),'.-','LineWidth',3,'MarkerSize',12)  
%     figure(105)
%     hold on
%     histogram(tetax(:,p), 24)
% %     histogram(tetay(:,p), 24)
% %     histogram(teta60(:,p), 24)
% end
% set(gca,'FontName','Cambria','FontSize',16);
% legend('T=1 K','T=5 K','T=10 K','T=20 K','T=30 K','T=35 K','T=50 K','T=60 K','T=75 K','T=100 K','T=150 K','T=200 K','T=250 K','T=300 K','T=400 K','T=500 K','Location','Best')
% xlabel(' cos( \theta ) ','Interpreter','tex' )
% ylabel('population','Interpreter','tex')