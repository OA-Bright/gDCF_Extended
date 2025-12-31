close all; clear all; clc

p = gcp('nocreate');
if isempty(p)
numworkers =8;
parpool('local', numworkers);
end

%Load in sample trajectory
%data size = [Nreadout, Nshots, Ntime_frame]
%trajectory is normalized within [-0.5 0.5]

%Some example caculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Example 1
data=load('./Radial_2D_36_spokes_GA.mat'); % enter filepath
k_space_traj = data.k_rad; 
size(k_space_traj)
%recon matrix size
N = 384;

% Delta_k approximated with Nyquist
    gDCF_1 = gDCF_extended_2D(k_space_traj,N,'Nyq',0,0, 'Slow'); %no nufft correction - slow calc
    gDCF_2 = gDCF_extended_2D(k_space_traj,N,'Nyq',0,0, 'Fast'); %no nufft correction - Fast calc
%Plot Traj and DCF for one timeframe 
figure(1)
plot(k_space_traj(:,:,1))
title(['k-space Traj'])

figure(2);
hold on
plot(gDCF_1(:,1,1),'r-','LineWidth',2.5,'DisplayName', 'gDCF-Nqy-Slow')
plot(gDCF_2(:,1,1),'b--','LineWidth',3.5,'DisplayName', 'gDCF-Nqy-Fast')
hold off
legend('show')
title(['gDCF-extended'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%Example 2
data=load('./Radial_2D_36_spokes_GA.mat'); % enter filepath
k_space_traj = data.k_rad; 
size(k_space_traj)
%recon matrix size
N = 384;


% Delta_k approximated with Nyquist
    gDCF_1 = gDCF_extended_2D(k_space_traj,N,'Nyq',0,0, 'Fast'); %no nufft correction
    gDCF_2 = gDCF_extended_2D(k_space_traj,N,'Nyq',0,1, 'Fast'); %nufft correction

 % Delta_k approximated with BL
 % Note BL type = 0 (full spokes) or 1 (center-out readout)
    gDCF_3 = gDCF_extended_2D(k_space_traj,N,'BL',0,0, 'Fast'); % no nufft correction; 
    gDCF_4 = gDCF_extended_2D(k_space_traj,N,'BL',0,1, 'Fast'); % nufft correction; 

%Plot Traj and DCF for one timeframe 
figure(3)
plot(k_space_traj(:,:,1))
title(['k-space Traj'])

figure(4);
hold on
plot(gDCF_1(:,1,1),'r-','LineWidth',2.5,'DisplayName', 'gDCF-Nqy')
plot(gDCF_2(:,1,1),'b--','LineWidth',3.5,'DisplayName', 'corr-gDCF-Nqy')
plot(gDCF_3(:,1,1),'g--','LineWidth',2.5,'DisplayName', 'gDCF-BL')
plot(gDCF_4(:,1,1),'m-','LineWidth',2,'DisplayName', 'corr-gDCF-BL')
hold off
legend('show')
title(['gDCF-extended'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Example 3
data=load('./Spiral_10Shots_18TF'); % enter filepath
k_space_traj = data.k_rad; 
size(k_space_traj)
%recon matrix size
N = 384;

% Delta_k approximated with Nyquist
    gDCF_1 = gDCF_extended_2D(k_space_traj,N,'Nyq',0,1, 'Fast'); %nufft correction

 % Delta_k approximated with BL
 % Note BL type = 0 (full spokes) or 1 (center-out readout)
    gDCF_2 = gDCF_extended_2D(k_space_traj,N,'BL',1,1, 'Fast'); % nufft correction; 

 % Delta_k approximated with Gen-FOV
    gDCF_3 = gDCF_extended_2D(k_space_traj,N,'GEN',0,1, 'Fast'); %nufft correction 

%Plot Traj and DCF for one timeframe 
figure(5)
plot(k_space_traj(:,:,1))
title(['k-space Traj'])

figure(6);
hold on
plot(gDCF_1(:,1,1),'r-','LineWidth',2.5,'DisplayName', 'corr-gDCF-Nqy')
plot(gDCF_2(:,1,1),'b--','LineWidth',3.5,'DisplayName', 'corr-gDCF-BL')
plot(gDCF_3(:,1,1),'g--','LineWidth',2.5,'DisplayName', 'corr-gDCF-GEN')
hold off
legend('show')
title(['gDCF-extended'])
