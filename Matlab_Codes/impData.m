% This code would import data from lammps log files into matlab 
clc, clear
delimiterIn = ' ';
headerlinesIn = 334 ; % 1-328 2-336 3-333 4-334 5-297 6-343
k=0;

% Data column order:
%1Step 2CPU 3PotEng 4KinEng 5Temp 6Lx 7Ly 8Press 
%9v_xc_x 10v_xc_y 11v_xc_z 12c_pe_c60 13c_lennard 14c_ke_c60 
%15v_vc_x 16v_vc_y 17v_vc_z
%18v_x1_x 19v_x1_y 20v_x1_z 21v_x2_x 22v_x2_y 23v_x2_z 24c_pe_sub 25c_ke_sub 
%26v_wc_x 27v_wc_y 28v_wc_z 29v_w12_x 30v_w12_y 31v_w12_z 32c_temp_c60 33c_temp_sub

%1Step 2KinEng 3PotEng 4Temp 5Press 6c_ke_nc 7c_lennard 8c_pe_nc 9v_xcm_x 10v_xcm_y 11v_xcm_z 
%12v_vcm_x v_vcm_y v_vcm_z v_wcm_x v_wcm_y v_wcm_z v_fcm_x v_fcm_y v_fcm_z 
%21c_msdd[1] c_msdd[2] c_msdd[3] v_x1_x v_x1_y v_x1_z v_x2_x v_x2_y v_x2_z 
%30v_wheel1_x v_wheel1_y v_wheel1_z v_wheel1_vx v_wheel1_vy v_wheel1_vz v_wheel1_ox v_wheel1_oy v_wheel1_oz 
%39v_wheel2_x v_wheel2_y v_wheel2_z v_wheel2_vx v_wheel2_vy v_wheel2_vz v_wheel2_ox v_wheel2_oy v_wheel2_oz 
%48v_wheel3_x v_wheel3_y v_wheel3_z v_wheel3_vx v_wheel3_vy v_wheel3_vz v_wheel3_ox v_wheel3_oy v_wheel3_oz 
%57v_wheel4_x v_wheel4_y v_wheel4_z v_wheel4_vx v_wheel4_vy v_wheel4_vz v_wheel4_ox v_wheel4_oy v_wheel4_oz 
%66v_chassi_x v_chassi_y v_chassi_z v_chassi_vx v_chassi_vy v_chassi_vz v_chassi_ox v_chassi_oy v_chassi_oz 
%75c_temp_nc c_temp_sub 

T=[5 10 30 50 75 100 150 200 300 400 500 600 700 800 900 1000]; 

for i=T         
    filename = sprintf('log.G12_1_T%d', i);
    eval(['temp = importdata(filename,delimiterIn,headerlinesIn)']);
    k=k+1
    imdata(:,:,k)=temp.data;
end
save allT1m.mat imdata



