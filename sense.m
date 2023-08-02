%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project       : Parallel Imaging - SENSE 
% Author        : Hakan Gulec
% Supervisor(s) : Berkin Bilgic & Yohan Jun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all

%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Description %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%

% This project implements SENSE(Sensitivity Encoding) parallel imaging
% technique for MRI image reconstruction. The data used in this code is
% provided as a .mat file. 

%-------------------------------------------------------------------------%


%% Load data


load sim_8ch_data.mat


%% Get the original image, coil sensitivity maps and coil images


% Ground truth image 
image_original = anatomy_orig / max(abs(anatomy_orig(:)));

% Coil sensitivity maps
coil_sensitivity_maps = b1 / max(abs(b1(:)));

% Number of coils, 8 for this data
number_of_coils = size(coil_sensitivity_maps, 3);

% MR images under each coil
image_coils = repmat(image_original, [1, 1, number_of_coils]) .* coil_sensitivity_maps;

% mosaic function is used to display 
mosaic(image_original       , 1, 1, 1, 'Original Image', [0  1], 0);
mosaic(coil_sensitivity_maps, 2, 4, 2, 'Coil Maps'     , [0 .5], 0);
mosaic(image_coils          , 2, 4, 3, 'Coil Images'   , [0 .5], 0);


%% Add Gaussian noise 


% Add noise to coil images 
noise_standard_deviation = 3e-3;       

% Noisy coil images
image_coils = image_coils + randn(size(image_coils)) * noise_standard_deviation;     


%% Create the aliased coil images


% Acceleration factors, I wrote this code only for R = 2, need some
% generalization
% R_x = 2;
R_y = 2;

% Create the aliased images when the acceleration factor R = 2
image_aliased = image_coils(1:(end / R_y), :, :) + image_coils((1 + end / R_y):end, :, :);

mosaic(image_aliased, 2, 4, 4, 'Aliased Images', [0 .5], 0);


%% Implement SENSE reconstruction


% Create the empty reconstructed image
reconstructed_image = zeros(size(image_original));

% Create the empty variables in the SENSE equation C_p * x = y where C_p is
% an N * R matrix containing the coil sensitivity values for this set of 
% pixels (N is the number of coils and R is the acceleration factor), y is 
% a vector of length N containing the aliased pixels for each coil, and x
% is the unknown but desired vector of R unaliased image pixels.
C_p = zeros(number_of_coils, R_y);
y   = zeros(number_of_coils, 1);

% Loop over the aliased pixels, for this data, it is 64*128
for i = 1:(size(image_original, 1) / R_y)

    for j = 1:size(image_original, 2)
          
        % Assign coil sensitivity values for the aliased pixels
        C_p(:, 1) = reshape(coil_sensitivity_maps(i          , j, :), 1, []).';
        C_p(:, 2) = reshape(coil_sensitivity_maps(i + end / 2, j, :), 1, []).';
        C_p       = abs(C_p)                                                  ;
        
        % Assign the image values for the aliased pixels
        y = reshape(image_aliased(i, j, :), 1, []).';
        y = abs(y)                                  ; 

        % Solve the equation x_hat = C_p_inverse*y
        temp = C_p \ y;
        temp = temp.' ;

        reconstructed_image(i          , j) = temp(1);
        reconstructed_image(i + end / 2, j) = temp(2);

    end

end


%% Display of the reconstructed image


mosaic(reconstructed_image, 1, 1, 5, 'Reconstructed Image', [0 1], 0);


%% Compute RMSE and display the reconstructed image


rmse_sense = 100 * norm(reconstructed_image(:) - image_original(:)) / norm(image_original(:));

mosaic(reconstructed_image, 1, 1, 6, ['SENSE rmse: ', num2str(rmse_sense)], [0, 1]);


