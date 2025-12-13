function [DCF_time] = gDCF_extended_2D(k,mat,Methods,BL_type,nufft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Written by Oluyemi Aboyewa
% % Version 2.0:  December, 2025.

%Input:
% %    - k space trajectory normalized between [-0.5,0.5]
% %    - k space dims = [Nreadout Nshots Ntimeframe]
% %    - mat = recon matrix size 
% %    - method = select which approximate dealta_k approach ["Nyquist", "Beruling-Landau", generalized-FOV"] 
% %    - BL_type = Readout type for BL ["center out", "full spoke"] 
% %    - nufft = if 1, gridding kernel corection is applied to gDCF

%Output:
% %    - DCF_time = gDCF  

% %    Ref: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
Ns = size(k, 1);
Nproj = size(k, 2);
DCF_time = zeros(size(k));

switch Methods
    case 'Nyq'
        UnitDistance = 1/ mat;
    case 'BL'
        if BL_type ==1 % Center out trajectory; e.g spiral
            lengths = linspace(0, mat, (mat/2) + 1);
            widths = linspace(0, mat, (mat/2) + 1);
        else        % full spokes  e.g radial projection, rosette
            lengths = linspace(0, mat, mat);
            widths = linspace(0, mat, mat);
        end
        [~, del_k] = LBD(k * mat, lengths, widths);
        UnitDistance = del_k/ mat;
    case 'GEN'
        UnitDistance = genFOV(k);

end

calc_prec =10^9;

input_X = real(k(:,:,1));
input_Y = imag(k(:,:,1));

[Ang, ~] = cart2pol(real(k), imag(k));
Ang = mod(Ang, 2 * pi);

DCF = zeros(Ns, Nproj);

parfor Selected_Proj = 1:Nproj
    xsel = input_X(:, Selected_Proj);
    ysel = input_Y(:, Selected_Proj);
    gDCF = zeros(Ns, 1);

    for Selected_Meas = 1:Ns

        Ref_X = xsel(Selected_Meas);
        Ref_Y = ysel(Selected_Meas);

        cond = input_X >= Ref_X - UnitDistance & input_X <= Ref_X + UnitDistance & ...
               input_Y >= Ref_Y - UnitDistance & input_Y <= Ref_Y + UnitDistance;
        indx_1 = find(cond);

        if ~isempty(indx_1)
            calc_dist = (UnitDistance - sqrt((Ref_X - input_X(indx_1)).^2 + (Ref_Y - input_Y(indx_1)).^2)) / UnitDistance;

            calc_dist(calc_dist < 0) = 0;
            gDCF(Selected_Meas) = 1 / sum(calc_dist);
        else
            gDCF(Selected_Meas) = 0;
        end
    end

    DCF(:, Selected_Proj) = round(gDCF*calc_prec)/calc_prec;
end

if nufft==1

    d = ones(size(input_X));
    b = ones(mat,mat);
    k_t = k(:,:,1);

    E = MCNUFFT_CPU_matrix(k_t,d,b);
    P = E.st{1,1}.p;

    DCF_vec = DCF(:);
    D_eff = P * (P' * DCF_vec);
    init_D =abs(D_eff);
    init_mean = mean(init_D);
    init_std =std(init_D);
    disp(['Initial Measure (Mean +/- std): ' num2str(init_mean) ' +/- ' num2str(init_std)]);

    W_corr = 1 ./ (abs(D_eff)+ eps);
    DCF_corrected_vec = DCF_vec .* W_corr;
    D_check = abs(P * (P' * DCF_corrected_vec));
    final_mean = mean(D_check);
    Density_std = std(D_check);
    disp(['Convergence Measure (Mean +/- std) should ~ 1: ' num2str(final_mean) ' +/- ' num2str(Density_std)]);

    DCF(:,:,1) = reshape(DCF_corrected_vec, size(k_t));
end

dcf_1 = abs(DCF(:,:,1) .* exp(1i * Ang(:,:,2:end)));
DCF_time(:,:,1) = DCF(:,:,1);
DCF_time(:,:,2:end) = dcf_1;
DCF_time( find( DCF_time == 0 ) ) = eps;
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);

end


function avg_dist = genFOV(k)
    % Input:
    % k: (Nsamples x Nproj), representing k-space points
    % Output:
    % avg_dist: Average of sorted distances to the six nearest neighbors
    
    kSpaceMatrix(:, :, 1) = real(k(:,:,1)); 
    kSpaceMatrix(:, :, 2) = imag(k(:,:,1)); %
    
    [Nsamples, Nproj, ~] = size(kSpaceMatrix);
    avgDistances = zeros(Nsamples * Nproj, 1);
    
    parfor pointIndex = 1:(Nsamples * Nproj)
        %sample and projection indices
        sample1 = mod(pointIndex - 1, Nsamples) + 1;
        proj1 = ceil(pointIndex / Nsamples);
        point1 = squeeze(kSpaceMatrix(sample1, proj1, :));  % Current point
        distances = zeros(Nsamples * (Nproj - 1), 1);
        count = 0;
        
        % Compare with points in other projections
        for proj2 = 1:Nproj
            if proj1 ~= proj2
                points2 = squeeze(kSpaceMatrix(:, proj2, :));
                dists = sqrt(sum((points2 - point1').^2, 2));
                distances(count + 1:count + Nsamples) = dists;
                count = count + Nsamples;
            end
        end
        
        sortedDistances = sort(distances(1:count)); % Sort distances 
        avgDistances(pointIndex) = mean(sortedDistances(1:min(6, count))); %select 6 nearest neighbors
    end
    
    avgDistances = avgDistances(avgDistances > 0);
    sortedAvgDistances = sort(avgDistances);
    avg_dist = mean(sortedAvgDistances);
end

function [D_B, delt_k] = LBD(sample_points, lengths, widths)
    
    % lengths: Vector of length values corresponding to each rectangle
    % widths: Vector of width values corresponding to each rectangle

    sample_x =real(sample_points(:,:,1));
    sample_y =imag(sample_points(:,:,1));
    assert(length(lengths) == length(widths), 'Lengths and widths must have the same number of elements.');
    
    num_rectangles = length(lengths);
    densities = zeros(num_rectangles-1, 1);
    np= zeros(num_rectangles-1, 1);
    na = zeros(num_rectangles-1, 1);

    for i = 1:num_rectangles-1
        L = lengths(i+1);
        W = widths(i+1);
        
        % Count points in a rectangle centered at the origin
        indices  = (sample_x >= -L/2 & sample_x <= L/2) & (sample_y >= -W/2 & sample_y <= W/2);
        N_LW =sum(indices(:));
        
        A_LW = L*W; 
        np(i) = N_LW;
        na(i) =A_LW;
        
        densities(i) = N_LW / A_LW; % number of points/ area of rectangle
    end

    % Return the lower Beurling density (limit inferior)
    non_zeroDB =densities(densities>=1);
    D_B = min(non_zeroDB(:));
    delt_k = 1/sqrt(D_B);
end

function  [res] = MCNUFFT_CPU_matrix(k,w,b1)

% Multicoil NUFFT operator
% Based on the NUFFT toolbox from Jeff Fessler and the single-coil NUFFT
% operator from Miki Lustig
% Input
% k: k-space trajectory
% w: density compensation
% b1: coil sensitivity maps
%
% Li Feng & Ricardo Otazo, NYU, 2012
%Modified by Oluyemi Aboyewa, 2025
    
Nd = size(b1(:,:,1));
Jd = [6,6];
Kd = floor([Nd*1.5]);
n_shift = Nd/2;
for tt=1:size(k,3)
	kk=k(:,:,tt);
	om = [real(kk(:)), imag(kk(:))]*2*pi;
	res.st{tt} = nufft_init(om, Nd, Jd, Kd, n_shift,'kaiser');
end

end
