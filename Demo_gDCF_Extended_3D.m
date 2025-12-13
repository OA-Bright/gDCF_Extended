close all; clear all; clc

%Load in sample trajectory
%[Nsample,Nshot,Ndim,Nt] =size(k) %Ndim = 3 for x,y,z,Nt= timeframe
% k space trajectory normalized between [-0.5,0.5,0.5]

%Example 1
data=load('./Radial_3D_103Shots.mat'); % enter filepath
k_space_traj = data.proj_grid; 
size(k_space_traj)
%recon matrix size
N = 128;

%compute DCF
gDCF_1 = gDCF_extended_3D(k_space_traj,N,0); %no nufft corection
gDCF_2 = gDCF_extended_3D(k_space_traj,N,1); %nufft corection

%Plot DCF for one time frame
figure(1); plot3( k_space_traj(:,:,1), k_space_traj(:,:,2), k_space_traj(:,:,3)) 
title(['k-space Traj'])

figure(2);
hold on
plot(gDCF_1(:,1,1),'r-','LineWidth',2.5,'DisplayName', 'gDCF-Nqy')
plot(gDCF_2(:,1,1),'b--','LineWidth',2.5,'DisplayName', 'corr-gDCF-Nqy')
hold off
legend('show')
title(['gDCF-extended'])

%Example 2
data=load('./Cones_106Shots.mat'); % enter filepath
k_space_traj = data.proj_grid; 
size(k_space_traj)
%recon matrix size
N = 128;

%compute DCF
gDCF_1 = gDCF_extended_3D(k_space_traj,N,0); %no nufft corection
gDCF_2 = gDCF_extended_3D(k_space_traj,N,1); %nufft corection

%Plot DCF for one time frame
figure(3); plot3( k_space_traj(:,:,1), k_space_traj(:,:,2), k_space_traj(:,:,3)) 
title(['k-space Traj'])

figure(4);
hold on
plot(gDCF_1(:,1,1),'r-','LineWidth',2.5,'DisplayName', 'gDCF-Nqy')
plot(gDCF_2(:,1,1),'b--','LineWidth',2.5,'DisplayName', 'corr-gDCF-Nqy')
hold off
legend('show')
title(['gDCF-extended'])
