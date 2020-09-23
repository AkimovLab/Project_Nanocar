clc, clear
load ../allT1.mat
T=[5 10 30 50 75 100 150 200 300 400 500 600 700 800 900 1000];
%  1 2  3  4  5  6   7   8   9   10  11  12  13  14  15  16


% Data Coupling
% size1=size(imdata,1);
% size2=size(imdata,2);
% Ignore=2002;
% for i=1:size1
%     Data(i,:,:)=reshape(imdata(Ignore:end,:,i),size2,[]);
% end

% vars = {'Data1'};
% clear(vars{:})

%1Step 2KinEng 3PotEng 4Temp 5Press 6c_ke_nc 7c_lennard 8c_pe_nc 9v_xcm_x 10v_xcm_y 11v_xcm_z 
%12v_vcm_x v_vcm_y v_vcm_z v_wcm_x v_wcm_y v_wcm_z v_fcm_x v_fcm_y v_fcm_z 
%21c_msdd[1] c_msdd[2] c_msdd[3] v_x1_x v_x1_y v_x1_z v_x2_x v_x2_y v_x2_z 
%30v_wheel1_x v_wheel1_y v_wheel1_z v_wheel1_vx v_wheel1_vy v_wheel1_vz v_wheel1_ox v_wheel1_oy v_wheel1_oz 
%39v_wheel2_x v_wheel2_y v_wheel2_z v_wheel2_vx v_wheel2_vy v_wheel2_vz v_wheel2_ox v_wheel2_oy v_wheel2_oz 
%48v_wheel3_x v_wheel3_y v_wheel3_z v_wheel3_vx v_wheel3_vy v_wheel3_vz v_wheel3_ox v_wheel3_oy v_wheel3_oz 
%57v_wheel4_x v_wheel4_y v_wheel4_z v_wheel4_vx v_wheel4_vy v_wheel4_vz v_wheel4_ox v_wheel4_oy v_wheel4_oz 
%66v_chassi_x v_chassi_y v_chassi_z v_chassi_vx v_chassi_vy v_chassi_vz v_chassi_ox v_chassi_oy v_chassi_oz 
%75c_temp_nc c_temp_sub 

Data=imdata;

Index_v_x=[33 42 51 60];                                      % Index of wheel velocity_x
Index_v_y=[34 43 52 61];                                      % Index of wheel velocity_y
Index_v_z=[35 44 53 62];                                      % Index of wheel velocity_z

Index_w_x=[36 45 54 63];                                      % Index of wheel omega_x
Index_w_y=[37 46 55 64];                                      % Index of wheel omega_y
Index_w_z=[38 47 56 65];                                      % Index of wheel omega_z

nsub=length(T);
col=jet(nsub);

%% Correlation of wheel's rotaional speed

N=length(Index_w_x);
N2=sqrt(3*N*(N-1)/2);

for k=1:length(T)
    Correl=[];
    p_value=[];
    for i=1:N
        for j=1:N
            if j>i
                x_i=reshape(Data(:,Index_w_x(i),k),[],1);
                x_j=reshape(Data(:,Index_w_x(j),k),[],1);
                y_i=reshape(Data(:,Index_w_y(i),k),[],1);
                y_j=reshape(Data(:,Index_w_y(j),k),[],1);
                z_i=reshape(Data(:,Index_w_z(i),k),[],1);
                z_j=reshape(Data(:,Index_w_z(j),k),[],1);
                
                [r_x,p_x]=corr(x_i,x_j);
                [r_y,p_y]=corr(y_i,y_j);
                [r_z,p_z]=corr(z_i,z_j);
                Correl=[Correl r_x r_y r_z];
                p_value=[p_value p_x p_y p_z];
                
            end
        end
    end
%     Norm_Corr(k)=   norm(Correl)/N2;
%     Norm_P_Value(k)=norm(p_value)/N2;
    Mean_Corr(k)=mean(abs(Correl));

end

% display(vpa(Norm_Corr',4))
% display(vpa(Norm_P_Value',4))

figure(10)
plot(T,Mean_Corr,'.','LineWidth',2,'MarkerSize',25)
hold on
xlabel('Temperature (K)','Interpreter','latex')
ylabel(' Mean Correlation ','Interpreter','latex')
set(gca,'FontName','Cambria','FontSize',20);

% legend('NC/SLG','NC/DLG','NC/FLG','Location','Best')
% legend('NC/SLG','NC/DLG','NC/FLG','NT/SLG','NT/DLG','NT/FLG','Location','Best')

fprintf('Mean Correlation =')
fprintf('% 10.4g',Mean_Corr')
fprintf('\n\n')


%% Wheels Slip Ratio

N_wheel=size(Index_v_x,2);                                          % No. of wheels

for k=1:length(T)
    
    SlipRatio=[];
    
    for i=1:N_wheel
        
        v_x=reshape(Data(:,Index_v_x(i),k),[],1);
        v_y=reshape(Data(:,Index_v_y(i),k),[],1);
        v_z=reshape(Data(:,Index_v_z(i),k),[],1);
        
        w_x=reshape(Data(:,Index_w_x(i),k),[],1);
        w_y=reshape(Data(:,Index_w_y(i),k),[],1);
        w_z=reshape(Data(:,Index_w_z(i),k),[],1);
        
        V=[v_x v_y v_z];
        W=[w_x w_y w_z];
        
        R=zeros(size(W));
        R(:,3)=3.5;
        WxR=cross(W,R);
        
        V_xy  =  V(:,1:2);
        WxR_xy=WxR(:,1:2);
        
        mean_speed=sqrt(mean(sum(V_xy.^2,2)));
        Index_Non_zero=find( sqrt(sum(V_xy.^2,2)) > mean_speed/10);
        WxR_xy_NZ=WxR_xy(Index_Non_zero,:);
        V_xy_NZ  =  V_xy(Index_Non_zero,:);
        
        slipratio=sqrt( sum((WxR_xy_NZ-V_xy_NZ).^2,2) ./ sum(V_xy_NZ.^2,2) );
        SlipRatio=[SlipRatio ; slipratio];
    end
    
    RMS_Slip_Ratio(k)=rms(SlipRatio);
    
end

fprintf('Slip Ratio=\n')
fprintf('% 5.4f   ',RMS_Slip_Ratio)
fprintf('\n')

figure(1);
plot(T,RMS_Slip_Ratio,'LineWidth',1.5)
xlabel('Temperature (K)','FontName','Cambria','FontSize',10)
ylabel('Slip Ratio','FontName','Cambria','FontSize',10)
grid on
hold on
set(gca,'FontName','Cambria','FontSize',16);

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

%% Rotation diffusion coefficient (Phys. Chem. Chem. Phys., 2013, 15, 845--849)

% figure(10)
Nseg=60;                                       % Segments No.
Nt=length(TimeI);                               % No. total data
Lseg=floor(Nt/Nseg);                            % Length of segment
TimeSeg=TimeI(1:Lseg);

for k=1:length(T)
    
    msd=zeros(Lseg,1);                          % Mean of MSD of all wheels @ temp k
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
        msd=msd+msd_temp;
    end
    msd=msd/N_wheel;
    
    ft = fittype('Poly1'); opts = fitoptions( ft ); opts.Lower = [0 0]; opts.Upper = [1000 0];
    fitresult=fit(TimeSeg, msd, ft, opts ); p1=fitresult.p1;    p2=fitresult.p2;
    
    MSD(k,:)=msd;
    P(k,:)=[p1 p2];
    
%     subplot(4,4,k)
%     plot(TimeSeg,msd,'LineWidth',2,'DisplayName','sys5')
%     hold on
% %     plot(TimeSeg,p1*TimeSeg+p2,'Color',[1 0 0],'LineWidth',1,'DisplayName','Fitted line')
%     title (['T=' num2str(T(k)) 'K'])
%     xlabel('Time (ps)')
%     ylabel('MSD (rad^2)')
%     legend('show','Location','Best');
%     set(gca,'FontName','Cambria','FontSize',12);
    
end


D=P(:,1)/4;
display(vpa(D,4))

figure(2)
plot(T,D,'LineWidth',2)
hold on
xlabel('Temperature (K)','FontName','Cambria','FontSize',18)
ylabel('D_{rot} (rad^2/ps)','FontName','Cambria','FontSize',18)
set(gca,'FontName','Cambria','FontSize',18);


%% Arrhenious
figure(3)
T_1= 1./T(2:end)';
lnDm= log(abs(D(2:end)));
plot(T_1,lnDm,'o','LineWidth',2,'MarkerSize',6)
hold on
xlabel('$ 1/T (K^{-1}) $','Interpreter','latex')
ylabel('ln(D_{rot}) (rad^2/ps)')
set(gca,'FontName','Cambria','FontSize',20);

% set(gca,'FontName','Cambria','FontSize',18,'YTick',[0.000001 0.00001 0.0001 0.001 0.01 0.1 1 10],'YScale','log','YMinorTick','on');
% ylim([0.0000001 0.01])

b1=polyfit(T_1(1:5),lnDm(1:5),1) %%
yhat1=polyval(b1,T_1(1:7));
plot(T_1(1:7),yhat1,'k--','MarkerSize',55);

%% Ea

F=abs(D);

kB=8.62*10^-2; %meV

LG=8; %1 4 7  2 5 8
figure(20)
bar(LG,-kB*b1(1),'r')
hold on
% errorbar(LG,-kB*b1(1),-kB*slopestd,'gx')
set(gca,'FontName','Cambria','FontSize',20);





