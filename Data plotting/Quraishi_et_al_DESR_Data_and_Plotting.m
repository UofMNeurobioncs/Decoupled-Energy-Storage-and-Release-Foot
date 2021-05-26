%% Section 0: Cleaning work environment

clc		% clears command window
clear all	% clears workspace
close all	% closes all open figures


%% Section 1: Theoretical desired torque-angle (TA) curves

angle_res = 0.005;				% Resolution, step size between angle data points
theta_total = pi/180*(-40:angle_res:40);	% Total angle range of motion for theoretical TA curve, in radians.


%% Section 2: Actual TA values measured on the Dynamometer

% Note that the dorsal and plantar curves are combined in a single array

load('DESR_Characterization_Data.mat')	% Loading the Dynamometer results           
[b,a] = butter(3,15/500);		% Butterworth filter settings

dynamometer_stiffness = 218*180/pi;					% Stiffness of the dynamometer, known prior to measurement
angle_dynamometer = DESR_N.x33.angleJIM;				% Ankle angle data [rad] from dynamometer at position 33 mm
torque_dynamometer = DESR_N.x33.torque;					% Torque data [N] at position 33mm, caused by frame and spring
torque_dynamometer_filtered = filtfilt(b,a,torque_dynamometer);		% Filtered torque data
torque_dynamometer_corrected = torque_dynamometer_filtered + 1.2;	% Corrected offset torque data, manual increase in torque, so the resting position is at 0N

angle_dynamometer_corrected = angle_dynamometer - torque_dynamometer_corrected/dynamometer_stiffness;	% Corrected angle rotation of the prosthesis ankle, accounted for elasticity of dynamometer


%% Section 3: Plotting desired TA curves and measured results

% In the 'DESR_Characterization_Data.dat' dataset, the desired torque-angle 
% curves of the 'plantar' and 'dorsal' cam profile are included for
% reference as 'torque_desired_plantar' and 'torque_desired_dorsal'

% New figure to plot the desired and measured TA curves
figure     
hold on

plot(theta_total/pi*180,torque_desired_plantar,'Linewidth',6,'Color', [88 164 176]/255,'LineStyle','--') 
plot(theta_total/pi*180,torque_desired_dorsal,'Linewidth',6,'Color', [255 136 0]/255,'LineStyle','--') 
plot(angle_dynamometer_corrected*(180/pi),torque_dynamometer_corrected,'linewidth',6 ,'Color',[55 94 151]/255) 

set(gcf,'color','w');
set(gca,'FontSize',25)
set(gca,'linewidth',2)
xlabel('Ankle Angle [deg]')
ylabel('Ankle Torque [N.m]')
xlim([-8 15])
ylim([-60 120])
legend('Desired plantar curve', 'Desired dorsal curve','Actual measurement','FontSize',20,'Location','northwest')
title('Two cam profiles')


%% Section 4: Energy calculations for theoretical TA curves

max_dorsiflexion = 10;	% Degrees
max_plantarflexion = 5;	% Degrees

[~,index_minimum_plantar_torque] = min(abs(torque_desired_plantar));	% Find the array index for which the torque is closest to zero
[~,index_minimum_dorsal_torque] = min(abs(torque_desired_dorsal));	% Find the array index for which the torque is closest to zero
[~,index_maximum_dorsal_torque] = min(abs(theta_total-max_dorsiflexion/180*pi));	% Find array index at which the angle is equal to max_dorsiflexion

% Energy values in joule
energy_theoretical_controlled_dorsiflexion = max(abs(cumtrapz(theta_total(index_minimum_dorsal_torque:index_maximum_dorsal_torque),torque_desired_dorsal(index_minimum_dorsal_torque:index_maximum_dorsal_torque))));
energy_theoretical_pushoff = max(abs(cumtrapz(theta_total(index_minimum_plantar_torque:index_maximum_dorsal_torque),torque_desired_plantar(index_minimum_plantar_torque:index_maximum_dorsal_torque))));
energy_theoretical_midstance = max(abs(cumtrapz(theta_total(1:index_minimum_plantar_torque),torque_desired_plantar(1:index_minimum_plantar_torque)))) - max(abs(cumtrapz(theta_total(1:index_minimum_dorsal_torque),torque_desired_dorsal(1:index_minimum_dorsal_torque)))); 

% Energy value in percentage [%]
energy_theoretical_recycled = (energy_theoretical_pushoff-energy_theoretical_controlled_dorsiflexion)/energy_theoretical_controlled_dorsiflexion*100; 


%% Section 5: Energy calculations for measured dual cam profiles

% Largest negative and positive angles
angle_plantar_min_index = find(angle_dynamometer_corrected == min(angle_dynamometer_corrected));
angle_plantar_max_index = find(angle_dynamometer_corrected == max(angle_dynamometer_corrected));

% Points where the curve intersects zero torque (x-axis)
[~,index_minimum_plantar_torque]= min(abs(torque_dynamometer_corrected(angle_plantar_min_index+1:angle_plantar_max_index)-0));
[~,index_maximum_dorsal_torque]= min(abs(torque_dynamometer_corrected(angle_plantar_max_index+1:end)-0));
   
% Array containing energy values at all measured points
energy_dualcam_controlled_dorsiflexion = max(abs(cumtrapz(angle_dynamometer_corrected(1:angle_plantar_min_index),torque_dynamometer_corrected(1:angle_plantar_min_index)))) + max(abs(cumtrapz(angle_dynamometer_corrected(angle_plantar_max_index+index_maximum_dorsal_torque:end),torque_dynamometer_corrected(angle_plantar_max_index+index_maximum_dorsal_torque:end))));
energy_dualcam_pushoff = max(abs(cumtrapz(angle_dynamometer_corrected(angle_plantar_min_index+1:angle_plantar_min_index+index_minimum_plantar_torque),torque_dynamometer_corrected(angle_plantar_min_index+1:angle_plantar_min_index+index_minimum_plantar_torque))));
energy_dualcam_midstance = energy_dualcam_controlled_dorsiflexion - energy_dualcam_pushoff;

% Max energy values: stored, released and recycled
energy_max_dualcam_controlled_dorsiflexion = max(abs(cumtrapz(angle_dynamometer_corrected(angle_plantar_min_index+index_minimum_plantar_torque+1:angle_plantar_max_index),torque_dynamometer_corrected(angle_plantar_min_index+index_minimum_plantar_torque+1:angle_plantar_max_index))));
energy_max_dualcam_pushoff = max(abs(cumtrapz(angle_dynamometer_corrected(angle_plantar_max_index+1:angle_plantar_max_index+index_maximum_dorsal_torque),torque_dynamometer_corrected(angle_plantar_max_index+1:angle_plantar_max_index+index_maximum_dorsal_torque))));
energy_max_dualcam_midstance = (energy_max_dualcam_pushoff-energy_max_dualcam_controlled_dorsiflexion)/energy_max_dualcam_controlled_dorsiflexion*100;

%% Section 6: Energy calculations for measured plantar cam profile (without switching)

angle_single_plantar = DESR_N_CAM1.x33.angleJIM;
torque_single_plantar = DESR_N_CAM1.x33.torque; 
torque_single_plantar_filtered = filtfilt(b,a,torque_single_plantar); 
torque_plantar_corrected = torque_single_plantar_filtered + 1.3;

% stiffness of the dynamometer doesn't change, see section 2 for the value
angle_plantar_corrected = angle_single_plantar - torque_plantar_corrected/dynamometer_stiffness;	%angle_single1 rotation of the VSPA, corrected for elasticity of dynamometer

figure
plot(angle_single_plantar, torque_plantar_corrected, 'linewidth',2)
grid on
xlabel('Angle (rad)')
ylabel('Torque (Nm)')

% Largest negative and positive angles 
angle_plantar_min_index = find(angle_plantar_corrected==min(angle_plantar_corrected));
angle_plantar_max_index = find(angle_plantar_corrected==max(angle_plantar_corrected));

% Points where the curve intersects zero torque (x-axis)
[~,index_minimum_plantar_torque]=min(abs(torque_plantar_corrected(angle_plantar_min_index+1:angle_plantar_max_index)-0));
[~,index_maximum_dorsal_torque]=min(abs(torque_plantar_corrected(angle_plantar_max_index+1:end)-0));
    
% Array containing energy values at all measured points
energy_plantarcam_controlled_dorsiflexion = max(abs(cumtrapz(angle_plantar_corrected(1:angle_plantar_min_index),torque_plantar_corrected(1:angle_plantar_min_index)))) + max(abs(cumtrapz(angle_plantar_corrected(angle_plantar_max_index+index_maximum_dorsal_torque:end),torque_plantar_corrected(angle_plantar_max_index+index_maximum_dorsal_torque:end))));
energy_plantarcam_pushoff = max(abs(cumtrapz(angle_plantar_corrected(angle_plantar_min_index+1:angle_plantar_min_index+index_minimum_plantar_torque),torque_plantar_corrected(angle_plantar_min_index+1:angle_plantar_min_index+index_minimum_plantar_torque))));
energy_plantarcam_midstance = energy_plantarcam_controlled_dorsiflexion - energy_plantarcam_pushoff;

% Max energy values: stored, released and recycled
energy_max_plantarcam_controlled_dorsiflexion = max(abs(cumtrapz(angle_plantar_corrected(angle_plantar_min_index+index_minimum_plantar_torque+1:angle_plantar_max_index),torque_plantar_corrected(angle_plantar_min_index+index_minimum_plantar_torque+1:angle_plantar_max_index))));
energy_max_plantarcam_pushoff = max(abs(cumtrapz(angle_plantar_corrected(angle_plantar_max_index+1:angle_plantar_max_index+index_maximum_dorsal_torque),torque_plantar_corrected(angle_plantar_max_index+1:angle_plantar_max_index+index_maximum_dorsal_torque))));
energy_max_plantarcam_midstance = (energy_max_plantarcam_pushoff-energy_max_plantarcam_controlled_dorsiflexion)/energy_max_plantarcam_controlled_dorsiflexion*100;

%% Section 7: Energy calculations for measured dorsal cam profile (without switching)

angle_single_dorsal = DESR_N_CAM2.x33.angleJIM;
torque_single_dorsal = DESR_N_CAM2.x33.torque; 
torque_single_dorsal_filtered = filtfilt(b,a,torque_single_dorsal);
torque_single_dorsal_corrected = torque_single_dorsal_filtered -1;

% stiffness of the dynamometer doesn't change, see section 2 for the value
angle_single_dorsal_corrected = angle_single_dorsal - torque_single_dorsal_corrected/dynamometer_stiffness;

figure
plot(angle_single_dorsal, torque_single_dorsal_corrected, 'linewidth', 2)
grid on
xlabel('Angle (rad)')
ylabel('Torque (Nm)')

% Largest negative and positive angles 
angle_plantar_min_index = find(angle_single_dorsal_corrected==min(angle_single_dorsal_corrected));
angle_plantar_max_index = find(angle_single_dorsal_corrected==max(angle_single_dorsal_corrected));

% Points where the curve intersects zero torque (x-axis)
[~,index_minimum_plantar_torque]=min(abs(torque_single_dorsal_corrected(angle_plantar_min_index+1:angle_plantar_max_index)-0));
[~,index_maximum_dorsal_torque]=min(abs(torque_single_dorsal_corrected(angle_plantar_max_index+1:end)-0));

% Array containing energy values at all measured points
energy_dorsalcam_controlled_dorsiflexion = max(abs(cumtrapz(angle_single_dorsal_corrected(1:angle_plantar_min_index),torque_single_dorsal_corrected(1:angle_plantar_min_index)))) + max(abs(cumtrapz(angle_single_dorsal_corrected(angle_plantar_max_index+index_maximum_dorsal_torque:end),torque_single_dorsal_corrected(angle_plantar_max_index+index_maximum_dorsal_torque:end))));
energy_dorsalcam_pushoff = max(abs(cumtrapz(angle_single_dorsal_corrected(angle_plantar_min_index+1:angle_plantar_min_index+index_minimum_plantar_torque),torque_single_dorsal_corrected(angle_plantar_min_index+1:angle_plantar_min_index+index_minimum_plantar_torque))));
energy_dorsalcam_midstance = energy_dorsalcam_controlled_dorsiflexion - energy_dorsalcam_pushoff;

% Max energy values: stored, released and recycled
energy_max_dorsalcam_controlled_dorsiflexion = max(abs(cumtrapz(angle_single_dorsal_corrected(angle_plantar_min_index+index_minimum_plantar_torque+1:angle_plantar_max_index),torque_single_dorsal_corrected(angle_plantar_min_index+index_minimum_plantar_torque+1:angle_plantar_max_index))));
energy_max_dorsalcam_pushoff = max(abs(cumtrapz(angle_single_dorsal_corrected(angle_plantar_max_index+1:angle_plantar_max_index+index_maximum_dorsal_torque),torque_single_dorsal_corrected(angle_plantar_max_index+1:angle_plantar_max_index+index_maximum_dorsal_torque))));
energy_max_dorsalcam_midstance = (energy_max_dorsalcam_pushoff-energy_max_dorsalcam_controlled_dorsiflexion)/energy_max_dorsalcam_controlled_dorsiflexion*100;

%% Section 8: Table generation of relevant results


% Column with descriptions
Energy = {'Energy stored at midstance (J)', 'Energy stored over dorsiflexion (J)', 'Energy returned during pushoff (J)', 'Energy returned (%)'}';

% Columns with values, theoretical and/or measured
DualCam_Theoretical = [energy_theoretical_midstance, energy_theoretical_controlled_dorsiflexion, energy_theoretical_pushoff, energy_theoretical_recycled]';
DualCam_Measured = [energy_dualcam_midstance, energy_max_dualcam_controlled_dorsiflexion, energy_max_dualcam_pushoff, energy_max_dualcam_midstance]'; 
PlantarflexionCam_Measured = [nan, energy_max_plantarcam_controlled_dorsiflexion, energy_max_plantarcam_pushoff, energy_max_plantarcam_midstance]';
DorsiflexionCam_Measured = [nan, energy_max_dorsalcam_controlled_dorsiflexion, energy_max_dorsalcam_pushoff, energy_max_dorsalcam_midstance]';

% Column classification: dual cam, single plantar cam or single dorsal cam
DualCam = table(DualCam_Theoretical,DualCam_Measured);
PlantarflexionCam = table(PlantarflexionCam_Measured);
DorsiflexionCam = table(DorsiflexionCam_Measured);

% Table with the desired results values.
% All rows are in Joule, except the last row, that is in percentage 
% The semicolon is removed, to print the result in the Command Window
ResultTable = table(Energy, DualCam, PlantarflexionCam, DorsiflexionCam) 
