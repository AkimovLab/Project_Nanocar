% pivoting rotation diffusion
clc, clear
load allT1.mat 
T=[5 10 30 50 75 100 150 200 300 400 500 600 700 800 900 1000];
%  1 2  3  4  5  6   7   8   9   10  11  12  13  14  15  16


s = 20000; % start point 10000

omega_z(:,:)=imdata(s:end,17,:);  % 72 73 74 for chassis

Nt = length(omega_z);             % number of data points
t = imdata(1:Nt,1)/1000;          % ps, total time
thermo=200;
dt = 1e-3;                        % 0.001 ps time between trajectory points


% integration 

rot_z=zeros(length(t),1);
Rot_z=zeros(length(t),length(T));

for j=1:length(T)

    rot_z(1)=0;
    for i=2:length(t)
        rot_z(i)=rot_z(i-1)+(t(i)-t(i-1))*(omega_z(i-1,j)+omega_z(i,j))/2;
    end
    Rot_z(:,j)=rot_z;
end

% TimeI=0:0.1:(length(imdata(:,1,1))-1)*0.1;        % ps
% TimeI=2*TimeI';
% 
% rot_z=zeros(length(TimeI),1);
% Rot_z=zeros(length(TimeI),length(T));
% 
% for k=1:length(T)
%     omega_z=imdata(:,17,k);    % 72 73 74 for chassis
%     rot_z(1)=0;
%     for j=2:length(TimeI)
%         rot_z(j)=rot_z(j-1)+(TimeI(j)-TimeI(j-1))*(omega_z(j-1)+omega_z(j))/2;
%     end
%     Rot_z(:,k)=rot_z;
% end


%% msd

Npar=60;                     % datapoints in each section  %1500 % 5000 so np =40
Nseg=floor(Nt/Npar);         % numver of particles   %33 % Length of segment
nt=1:Nseg;
time=dt*thermo*nt;              % time in (ps)

nn1=5;
nn2=Npar/nn1;

TimeSeg=t(1:Nseg);

R=zeros(Nseg,Npar);
msdk = zeros(Nseg,nn1,length(T)); msd = zeros(Nseg,length(T));
for j=1:length(T)
    for ts = 1:Nseg;
        for k=1:nn1
            for i = (k-1)*nn2:k*nn2-1
                msdk(ts,k,j)= msdk(ts,k,j)+((Rot_z(1+ts+i*Nseg,j)-Rot_z(1+i*Nseg,j))^2)/nn2;
                %
            end
            msd(ts,j)=mean(msdk(ts,:,j));
        end
    end
end

start=1000;

ss=time(start:end)-time(start);  %[1:Nseg];
for j=1:length(T)
    for k=1:nn1
        mm=msdk(start:end,k,j)-msdk(start,k,j);
        [m]=polyfit(ss',mm,1);
        ft = fittype( 'a*x', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( ft );
        opts.StartPoint = 450; %0.75
        [fitresult, gof] = fit( ss',mm, ft, opts );
        F(k,j)=fitresult.a/2;
        D(k,j)=m(1)/2;
    end
    Fmean(j)=mean(F(:,j));
    Fstd(j)=std(F(:,j));
    Dmean(j)=mean(D(:,j));
    Dstd(j)=std(D(:,j));
end

%% ploting MSD Curves
nsub=length(T);
col=jet(nsub); % parula, hsv, hot, pink, flag
               %  autumn, bone, colorcube, cool, copper, gray, 
               % jet, lines, prism, spring, summer, winter 

figure(12)
for p=1:6
    txt = ['T = ',num2str(T(p)), 'K'];
plot(time,msd(:,p),'.-','LineWidth',3,'MarkerSize',12,'color',col(p,:),'DisplayName',txt)
hold on
end
xlabel('t (ps)','Interpreter','latex')
ylabel('MSD ($rad^{2}$)','Interpreter','latex')
legend show
set(gca,'FontName','Cambria','FontSize',20);


%% D-T

figure(13)
% plot(T,F,'LineWidth',2)
errorbar(T,Fmean,Fstd,'linewidth',2,'MarkerSize',12)
hold on
xlabel('Temperature (K)','fontsize',20,'FontName','Cambria')
ylabel('D_{rot} (rad^2/ps)','fontsize',20,'FontName','Cambria')
set(gca,'FontName','Cambria','FontSize',20);

% Arrhenius
figure(15)
T_1= 1./T(2:end);
lnDm= log(abs(Fmean(2:end)));

plot(T_1,lnDm,'o','LineWidth',2,'MarkerSize',6)
hold on
xlabel('$ 1/T (K^{-1}) $','Interpreter','latex')
ylabel('ln(D_{rot}) (rad^2/ps)')
set(gca,'FontName','Cambria','FontSize',20);

b1=polyfit(T_1(1:4),lnDm(1:4),1) %%
yhat1=polyval(b1,T_1(1:7));
plot(T_1(1:7),yhat1,'k--','MarkerSize',55);

%% errorbar

F=abs(F);

figure(155)
for k=1:nn1
    lnD= log(F(k,2:end));
    plot(T_1, lnD ,'o','LineWidth',2,'MarkerSize',12)  
    hold on
    b2=polyfit(T_1(1:4),lnD(1:4),1); %%
    yhat2=polyval(b2,T_1(1:7));
    plot(T_1(1:7),yhat2,'k--','MarkerSize',55);
    slope(k)= b2(1);
end
slopem =mean(slope)
slopestd=std(slope)
sloperng=range(slope)
kB=8.62*10^-2; %meV

LG=8; %1 4 7  2 5 8
figure(20)
bar(LG,-kB*b1(1))
hold on
errorbar(LG,-kB*b1(1),-kB*slopestd,'rx')
set(gca,'FontName','Cambria','FontSize',20);

% legend('NC/SLG','NC/DLG','NC/FLG','NT/SLG','NT/DLG','NT/FLG','Location','Best')

    
