function [DCF_time] = gDCF_extended_3D(k,mat,nufft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Written by Oluyemi Aboyewa
% % Version 2.0:  December, 2025.
%Input:
% %    [Nsample,Nshot,Ndim,Nt] =size(k) %Ndim = 3 for x,y,z,Nt= timeframe
% %    - k space trajectory normalized between [-0.5,0.5,0.5]
% %    - mat = recon matrix size 
% %    - nufft = if 1, gridding kernel corection is applied to gDCF
%Output:
% %    DCF_time =gDCF

% %    Ref: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
Ns = size(k, 1);
Nproj = size(k, 2);
Nt = size(k,4); 

for ii=1:Nt
  fprintf( 'TimeFrame = %3d\n', ii );

  UnitDistance = 1/(mat); 

  input_X = k(:,:,1,ii);
  input_Y = k(:,:,2,ii);
  input_Z = k(:,:,3,ii);

  DCF = zeros(Ns, Nproj);

  parfor Selected_Proj = 1:Nproj
    xsel = input_X(:, Selected_Proj);
    ysel = input_Y(:, Selected_Proj);
    zsel = input_Z(:, Selected_Proj);

    gDCF = zeros(Ns, 1);

    for Selected_Meas = 1:Ns
      Ref_X = xsel(Selected_Meas);
      Ref_Y = ysel(Selected_Meas);
      Ref_Z = zsel(Selected_Meas);

      tmp_a = (input_X >= (Ref_X-UnitDistance)) & (input_X <= (Ref_X+UnitDistance));
      tmp_b = (input_Y >= (Ref_Y-UnitDistance)) & (input_Y <= (Ref_Y+UnitDistance));
      tmp_c = (input_Z >= (Ref_Z-UnitDistance)) & (input_Z <= (Ref_Z+UnitDistance));
      mask = tmp_a & tmp_b & tmp_c;
      indx_1 = find( mask > 0 );

      if ~isempty(indx_1)

                dist_sq = (Ref_X-input_X(indx_1)).^2 + (Ref_Y-input_Y(indx_1)).^2 + (Ref_Z-input_Z(indx_1)).^2;
                dist = sqrt(dist_sq);

                calc_dist = (UnitDistance - dist) / UnitDistance; 

        calc_dist(calc_dist < 0) = 0;

        gDCF(Selected_Meas) = 1 / (sum(calc_dist) + eps);
      else
        gDCF(Selected_Meas) = 0;
      end
    end

    DCF(:, Selected_Proj) = gDCF;
  end

 if nufft==1
  d = ones(Ns,Nproj);
  b = ones(mat, mat, mat);
  k_t = k(:,:,:,ii);

  E = MCNUFFT_3D_matrix(k_t,d,b); 
  P = E.st.p;
  DCF_vec = DCF(:);

  D_eff = P * (P' * DCF_vec);
  init_D =abs(D_eff);
  init_mean = mean(init_D);
  init_std =std(init_D);
  disp(['Initial Measure (Mean +/- std): ' num2str(init_mean) ' +/- ' num2str(init_std)]);

  W_corr = 1 ./ (init_D + eps);
  DCF_corrected_vec = DCF_vec .* W_corr;

  D_check = abs(P * (P' * DCF_corrected_vec));

  final_mean = mean(D_check);
  Density_std = std(D_check);

  disp(['Convergence Measure (Mean +/- std) should ~ 1: ' num2str(final_mean) ' +/- ' num2str(Density_std)]);

  DCF = reshape(DCF_corrected_vec, [Ns, Nproj]); 
end

DCF_time(:,:,:,ii) = DCF; 
end
DCF_time( find( DCF_time == 0 ) ) = eps;
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);

end

function  [res] = MCNUFFT_3D_matrix( traj, w, b1)
% Multicoil NUFFT operator
% Based on the NUFFT toolbox from Jeff Fessler and the single-coil NUFFT
% operator from Miki Lustig
% Input
% traj: k-space trajectory
% w: density compensation
% b1: coil sensitivity maps
%
% Li Feng & Ricardo Otazo, NYU, 2012
% Modified by Oluyemi Aboyewa, 2025

traj_x = col( traj(:,:,1) );
traj_y = col( traj(:,:,2) );
traj_z = col( traj(:,:,3) );

Nd = size(b1(:,:,:,1));
Jd = [6,6,6];
Kd = floor([Nd*1.5]);
n_shift = Nd/2;
om = [traj_x, traj_y, traj_z]*2*pi;
res.st = nufft_init(om, Nd, Jd, Kd, n_shift,'kaiser');

end
