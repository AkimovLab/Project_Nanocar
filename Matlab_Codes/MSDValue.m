% This code calculate the MSD value for each temperature by ensemble average 
clc, clear
load allT1.mat 
T=[5 10 30 50 75 100 150 200 300 400 500 600 700 800 900 1000];
%  1 2  3  4  5  6   7   8   9   10  11  12  13  14  15  16


s = 20000; % start point 10000

x(:,:)=imdata(s:end,9,:);    % throw off 5% of initial data points
y(:,:)=imdata(s:end,10,:);
Nt = length(x);              % number of data points
% t = imdata(1:Nt,1)/.1000;    % ps, total time
thermo=200;
dt = 1e-3;                   % 0.001 ps time between trajectory points
Tc=imdata(:,75,:);

Npar=60;                     % datapoints in each section  %1500 % 5000 so np =40
Nseg=floor(Nt/Npar);         % numver of particles   %33
%NT = floor(nData); %# for MSD, dt should be up to 1/4 of number of data points
nt=1:Nseg;
time=dt*thermo*nt;              % time in (ps)

nn1=5;
nn2=Npar/nn1;

msdk = zeros(Nseg,nn1,length(T)); msd = zeros(Nseg,length(T));
% X=zeros(Npar,Nseg); Y=zeros(Npar,Nseg);
start=1000;

%% msd

for j=1:length(T)               % calculate msd for all T's
    for ts = 1:Nseg;
        for k=1:nn1
            for i = (k-1)*nn2:k*nn2-1
                msdk(ts,k,j)= msdk(ts,k,j)+((x(1+ts+i*Nseg,j)-x(1+i*Nseg,j))^2+(y(1+ts+i*Nseg,j)-y(1+i*Nseg,j))^2)/nn2;
%               msdk2(ts,k,j)=((X(i+1,ts)-X(i+1,1))^2+(Y(i+1,ts)-Y(i+1,1))^2)/nn2+msdk2(ts,k,j);
            end
        end
        msd(ts,j)=mean(msdk(ts,:,j));
    end  
end



ss=time(start:end)-time(start);  %[1:Nseg];
for j=1:length(T)
    for k=1:nn1
        mm=msdk(start:end,k,j)-msdk(start,k,j);
        [m]=polyfit(ss',mm,1);
        ft = fittype( 'a*x', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( ft );
        opts.StartPoint = 450; %0.75
        [fitresult, gof] = fit( ss',mm, ft, opts );
        F(k,j)=fitresult.a/4;
        D(k,j)=m(1)/4;
    end
    Fmean(j)=mean(F(:,j));
    Fstd(j)=std(F(:,j));
    Dmean(j)=mean(D(:,j));
    Dstd(j)=std(D(:,j));
end



%% ploting MSD Curves
% % nsub=length(T);
% % col=jet(nsub); % parula, hsv, hot, pink, flag
% %                %  autumn, bone, colorcube, cool, copper, gray, 
% %                % jet, lines, prism, spring, summer, winter    
% % iall=1:16;
% % i1=1:8;
% % i2=8:18;
% % 
% % 
% % figure(11)
% % for p=iall
% %     txt = ['T = ',num2str(T(p)), 'K']
% %     plot(time,msd(:,p),'.-','LineWidth',3,'MarkerSize',12,'color',col(p,:),'DisplayName',txt)
% %     %M(p,1) = time(:)\msd(:,p);     % linear fitting with zero intercept
% %     % y_est = tt*M(p);
% %     % plot(tt, y_est, '-r')
% %     Tnow=round(imdata(end,75:76,p))
% %     hold on
% % end
% % xlabel('t (ps)','Interpreter','latex')
% % ylabel('MSD ($\AA^{2}$)','Interpreter','latex')
% % set(gca,'FontName','Cambria','FontSize',20);
% % legend show
% % 
% % figure(12)
% % for p=1:6
% %     txt = ['T = ',num2str(T(p)), 'K'];
% % plot(time,msd(:,p),'.-','LineWidth',3,'MarkerSize',12,'color',col(p,:),'DisplayName',txt)
% % hold on
% % end
% % xlabel('t (ps)','Interpreter','latex')
% % ylabel('MSD ($\AA^{2}$)','Interpreter','latex')
% % legend show
% % set(gca,'FontName','Cambria','FontSize',20);
% % 
% % figure(13)
% % for p=6:14
% %     txt = ['T = ',num2str(T(p)), 'K'];
% % plot(time,msd(:,p),'.-','LineWidth',3,'MarkerSize',12,'color',col(p,:),'DisplayName',txt)
% % hold on
% % end
% % xlabel('t (ps)','Interpreter','latex')
% % ylabel('MSD ($\AA^{2}$)','Interpreter','latex')
% % legend show
% % set(gca,'FontName','Cambria','FontSize',20);
% % 
% % 
% % %% log log msd
% % a=1;
% % b=3166;
% % lmsd=log10(msd);
% % ltime=log10(time);
% % 
% % for p=[2:16]
% %     txt = ['T = ',num2str(T(p)), 'K'];
% %     
% % %     figure(16)
% % %     loglog(time(a:end),msd(a:end,p),'.-','color',[rand,rand,rand],'MarkerSize',25,'DisplayName',txt)
% % %     hold on
% % 
% %     figure(17)
% %     plot(ltime(a:end),lmsd(a:end,p),'.-','color',[rand,rand,rand],'MarkerSize',25,'DisplayName',txt)
% %     hold on
% %     
% %     b=polyfit(ltime(a:end),lmsd(a:end,p)',1);
% %     Blog(p,:)=b;
% % end
% % % D=(10.^(Blog(:,2)))/4;
% % 
% % xlabel('log(t)(ps) ','Interpreter','latex')
% % ylabel('log(MSD)($\AA^{2}$) ','Interpreter','latex')
% % legend show
% % set(gca,'FontName','Cambria','FontSize',14);
% % 
% % % figure(16)
% % % xlabel('t (ps)','Interpreter','latex')
% % % ylabel('MSD ($\AA^{2}$)','Interpreter','latex')
% % % legend show
% % % set(gca,'FontName','Cambria','FontSize',12);
% % 
% % 
% % for p=[9]
% %     txt = ['T = ',num2str(T(p)), 'K'];
% %     figure(18)
% %     loglog(time,msd(:,p),'.-','color',[rand,rand,rand],'MarkerSize',25,'DisplayName',txt)
% %     hold on
% % end
% % xlabel('t (ps)','Interpreter','latex')
% % ylabel('MSD ($\AA^{2}$)','Interpreter','latex')
% % legend show
% % set(gca,'FontName','Cambria','FontSize',16);
% % 
% % % line([1.5 3],[3 6],'Color','black','LineStyle','--')
% % % line([1.5 3],[0.5 2],'Color','black','LineStyle','--')
% % 
% % %% calculation of alfa power
% % 
% % alfa=smooth(Blog(:,1),'sgolay');
% % figure(3)
% % hold on
% % plot(T, alfa,'o','LineWidth',0.5,'MarkerSize',12)
% % xlabel('Temperature (K)','Interpreter','latex')
% % ylabel('$$\alpha$$','Interpreter','latex')
% % set(gca,'FontName','Cambria','FontSize',20);
% % % line([30 30],[0 2],'Color','black','LineStyle','--')
% % % line([150 150],[0 2],'Color','black','LineStyle','--')
% % % axis([5 500 0.5 2])
% % 
% % 
% % figure(1)
% % errorbar(T,Fmean,Fstd,'linewidth',2,'MarkerSize',12)
% % % errorbar(T,Dmean,Dstd,'linewidth',3,'MarkerSize',10)
% % % plot(T,Dmean,'-bs','linewidth',3,'MarkerSize',12)
% % hold on
% % grid off
% % xlabel('Temperature (K)','fontsize',20,'FontName','Cambria')
% % ylabel('D (Å^{2}/ps)','fontsize',20,'FontName','Cambria')
% % set(gca,'FontName','Cambria','FontSize',20);



%% Arrhenius
figure(2)
T_1= 1./T(2:end); %%
lnDm= log(abs(Fmean(2:end))); %%
% lnDstd= log(Fstd(2:end));
plot(T_1, lnDm,'o','LineWidth',2,'MarkerSize',12)
hold on
% yy = smooth(lnD,'sgolay');
% plot(T_1, yy,'.-','LineWidth',4,'MarkerSize',16)
xlabel('$ 1/T (K^{-1}) $','Interpreter','latex')
ylabel('lnD ($\AA^{2}/ps$)','Interpreter','latex')
set(gca,'FontName','Cambria','FontSize',20);

b1=polyfit(T_1(1:3),lnDm(1:3),1) %%
yhat1=polyval(b1,T_1(1:7));
plot(T_1(1:7),yhat1,'k--','MarkerSize',55);
%at=9;
%b2=polyfit(T_1(at:end),lnDm(at:end),1);
%yhat2=polyval(b2,T_1(at-3:end));
%plot(T_1(at-3:end),yhat2,'k--','MarkerSize',55)

%% errorbar

F=abs(F);

figure(155)
for k=1:nn1
    lnD= log(F(k,2:end));
    plot(T_1, lnD ,'o','LineWidth',2,'MarkerSize',12)
    hold on
    b2=polyfit(T_1(1:3),lnD(1:3),1); %%
    yhat2=polyval(b2,T_1(1:7));
    plot(T_1(1:7),yhat2,'k--','MarkerSize',55);
    slope(k)= b2(1);
end
slopem =mean(slope)
slopestd=std(slope)
sloperng=range(slope)
kB=8.62*10^-2; %meV

LG=2; %1 4 7  2 5 8
figure(20)
bar(LG,-kB*b1(1))
hold on
errorbar(LG,-kB*b1(1),-kB*slopestd,'rx')
set(gca,'FontName','Cambria','FontSize',20);

% bC60=-18.2133;
% bar(3,-kB*bC60)

% legend('NC/SLG','NC/DLG','NC/FLG','Location','Best')
% legend('NC/SLG','NC/DLG','NC/FLG','NT/SLG','NT/DLG','NT/FLG','Location','Best')



