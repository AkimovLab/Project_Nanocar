% pivoting rotation diffusion
clc, clear
load allT6.mat 
T=[5 10 30 50 75 100 150 200 300 400 500 600 700 800 900 1000];
%  1 2  3  4  5  6   7   8   9   10  11  12  13  14  15  16

Data=imdata;

Index_v_x=[33 42 51 60];                                      % Index of wheel velocity_x
Index_v_y=[34 43 52 61];                                      % Index of wheel velocity_y
Index_v_z=[35 44 53 62];                                      % Index of wheel velocity_z

Index_w_x=[36 45 54 63];                                      % Index of wheel omega_x
Index_w_y=[37 46 55 64];                                      % Index of wheel omega_y
Index_w_z=[38 47 56 65];                                      % Index of wheel omega_z

s = 20000; % start point 10000

% omega_z(:,:)=imdata(s:end,17,:);  % 72 73 74 for chassis

% Nt = length(omega_z);             % number of data points
% t = imdata(1:Nt,1)/1000;          % ps, total time
thermo=200;
dt = 1e-3;                        % 0.001 ps time between trajectory points


% integration 
%% Rotation by integration of angular velocity

Index=[36 37 45 46 54 55 63 64];            % Index for integration of wheel rotation around x and y-axis
N_wheel=length(Index)/2;                                % No of wheel/2,    x and y omega

TimeI=0:0.1:(length(Data(:,1,1))-1)*0.1;                % ps
TimeI=2*TimeI';

rot=zeros(length(TimeI),1);
Rot=zeros(length(TimeI),length(T),length(Index));	% Rot(:,Temp,Wheel)

for i=1:length(Index)
    for k=1:length(T)
        omega=reshape(Data(:,Index(i),k),[],1);
        rot(1)=0;
        for j=2:length(TimeI)
            rot(j)=rot(j-1)+(TimeI(j)-TimeI(j-1))*(omega(j-1)+omega(j))/2;
        end
        Rot(:,k,i)=rot;                         
    end
end

% there is alternative!

%% Rotation diffusion coefficient (Phys. Chem. Chem. Phys., 2013, 15, 845--849)

% figure(1)


Nseg=60;                                       % Segments No.
Nt=length(TimeI);                               % No. total data
Lseg=floor(Nt/Nseg);                            % Length of segment
TimeSeg=TimeI(1:Lseg);
nt=1:Lseg;
time=dt*thermo*nt;              % time in (ps)


msdk=zeros(Lseg,length(T),length(Index));        % Mean of MSD of all wheels @ temp k
    

for k=1:length(T)
    msdt=zeros(Lseg,1);  
    for i=1:length(Index)
        
        Rot_temp=Rot(:,k,i);
        Rot_Temp=zeros(Lseg,Nseg);
        for j=1:Nseg
            Rot_Temp(1:Lseg,j)=Rot_temp((j-1)*Lseg+1:j*Lseg,1);
        end
        
        msd_temp=zeros(Lseg,1);
        for j=1:Lseg
            msd_temp(j,1)=mean((Rot_Temp(j,:)-Rot_Temp(1,:)).^2);
        end
        msdt=msdt+msd_temp;
        msdt=msdt/N_wheel;
        msdk(:,k,i)=msdt;
    end
%     msdt=msdt/N_wheel;
%     msdk=msdk/N_wheel;
    msd(:,k)=msdt;
end

% msd=mean(msdk,3);
   
start=1000;

ss=time(start:end)-time(start);  %[1:Nseg];
for j=1:length(T)
    for k=1:length(Index)
        mm=msdk(start:end,j,k)-msdk(start,j,k);
        [m]=polyfit(ss',mm,1);
        ft = fittype( 'a*x', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( ft );
        opts.StartPoint = 450; %0.75
        [fitresult, gof] = fit( ss',mm, ft, opts );
        F(k,j)=3*fitresult.a/2;
%         D(k,j)=m(1)/1;
    end
    Fmean(j)=mean(F(:,j));
    Fstd(j)=std(F(:,j));
%     Dmean(j)=mean(D(:,j));
%     Dstd(j)=std(D(:,j));
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

b1=polyfit(T_1(2:5),lnDm(2:5),1) %%
yhat1=polyval(b1,T_1(1:7));
plot(T_1(1:7),yhat1,'k--','MarkerSize',55);

%% errorbar

F=abs(F);

figure(155)
for k=1:length(Index)
    lnD= log(F(k,2:end));
    plot(T_1, lnD ,'o','LineWidth',2,'MarkerSize',12)  
    hold on
    b2=polyfit(T_1(2:5),lnD(2:5),1); %% 2:5 (30 to 100)
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
errorbar(LG,-kB*b1(1),-kB*slopestd,'gx')
set(gca,'FontName','Cambria','FontSize',20);

% legend('NC/SLG','NC/DLG','NC/FLG','NT/SLG','NT/DLG','NT/FLG','Location','Best')
% legend('NT/SLG','NT/DLG','NT/FLG','Location','Best')
    
