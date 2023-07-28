function dictionary = generate_dictionary_deterministic()
% Dictionary Generation for Microstructure Informed Myelin Mapping (MIMM)

% Mert Sisman 7/27/2023

rng('default') % for reproducibility 


matrix_size = 512;
%Set volume susceptibilities.
chi_iso = -0.1e-6;    %% -0.1 ppm / isotropic susceptibility of myelin sheaths.
chi_ani = -0.1e-6;    %% -0.1 ppm / anisotropic susceptibility of myelin sheaths.
chi_iron = 0.3e-6;  %% 0.3 ppm  / isotropic iron susceptibility.


m = 1;

for density = linspace(1e-3,1,10) %% 2DFD

% Define input parameters for 2D circle distribution generation.
S.circSize = [  matrix_size/180 matrix_size/160 matrix_size/140 matrix_size/120 matrix_size/100 matrix_size/80 matrix_size/60 matrix_size/40 matrix_size/24 matrix_size/20]; 
S.nSizes = NaN;
S.frameSize = [matrix_size*3/4 matrix_size*3/4];
S.edgeType = 1;
S.supressWarning = true;
S.density = density;
S.drawFrame = false;
%Circle distribution generation.
[circData, circHandles, frame, S] = bubblebath(S);
circData = sortrows(circData,3);

% Create a 2D grid.
lim = (matrix_size - 1) / 2;
x = -lim:lim;
y = -lim:lim;
[X,Y] = meshgrid(x,y);
Y = flipud(Y);




% Make sure that circle distribution which represents a fiber bundle
% resides inside a circular regon.
for l = length(circData):-1:1
    if circData(l,1)^2 + circData(l,2)^2 >  (matrix_size*3/8+circData(l,3))^2 
        circData(l,:) = []; 
    end    
end

for g_ratio = linspace(0.5,1,11) %g-ratio (ri/ro).
for theta =  0:5:90              %Fiber orientation.
for EID = 0:0.2:1                %Extracellular iron density.

% Initilizate circles mask and hollow circles mask that represent fibers and myelin
% 2D.
circles = zeros(matrix_size);
hollow_circles = zeros(matrix_size);

% Initilizate circle center mask.
centers = zeros(matrix_size,matrix_size,length(S.circSize));

for i = 1:length(circData)
    % Create a mask for the inner and outer circles of a single element.
    one_circle_outer = ((X-circData(i,1)).^2 + (Y-circData(i,2)).^2) < circData(i,3).^2;
    one_circle_inner = ((X-circData(i,1)).^2 + (Y-circData(i,2)).^2) < (circData(i,3).*g_ratio).^2;
    
    % Add the new circle to the total mask;
    circles = circles + one_circle_outer;
    hollow_circles = hollow_circles + one_circle_outer - one_circle_inner;
    
    % Determine the center of each circle and place it in the discrete
    % center mask.
    for j = 1:length(S.circSize)  
        if circData(i,3) == S.circSize(j)
            centers(round(circData(i,1)+matrix_size/2),round(circData(i,2)+matrix_size/2),j) = 1;
        end
    end
end


%%

% Initilizate 2D myelin related field distribution.
fields_multi_circle = zeros(matrix_size,matrix_size,length(S.circSize));

% Compute the field perturbation due to each hollow circle and shift it
% according the center location.
for j = 1:length(S.circSize)  
    field = compute_field(matrix_size,S.circSize(j),g_ratio, theta, chi_iso, chi_ani);
    fields_multi_circle(:,:,j) = conv2(centers(:,:,j),field,'same');
end

% Sum the field pertubations of each individual hollow circle.
total_field_2D = sum(fields_multi_circle,3);




%%
% Extention to 3D.

% Replicate 2D field, hollow circle mask and circle mask to 3D.
total_field_3D = repmat(total_field_2D, 1, 1, matrix_size); 
myelin_mask_resamp = repmat(hollow_circles, 1, 1, matrix_size); 
fiber_mask_resamp = repmat(circles, 1, 1, matrix_size); 

% Resample the 3D field and the masks to the defined fiber orientation
% around -y axis.
direction = [0 -1 0];
total_field_3D_resamp = imrotate3(total_field_3D, theta, direction,'linear','crop','FillValues',0);
myelin_mask_resamp = imrotate3(myelin_mask_resamp, theta, direction,'nearest','crop','FillValues',0);
fiber_mask_resamp = imrotate3(fiber_mask_resamp, theta, direction,'nearest','crop','FillValues',0);


% Extract the central ROI that the mGRE signal will be calculated. Crop the
% outer regions.
myelin_field = total_field_3D_resamp(matrix_size/4+1:end-matrix_size/4,matrix_size/4+1:end-matrix_size/4,matrix_size/4+1:end-matrix_size/4);
myelin_mask_resamp = myelin_mask_resamp(matrix_size/4+1:end-matrix_size/4,matrix_size/4+1:end-matrix_size/4,matrix_size/4+1:end-matrix_size/4);
fiber_mask_resamp = fiber_mask_resamp(matrix_size/4+1:end-matrix_size/4,matrix_size/4+1:end-matrix_size/4,matrix_size/4+1:end-matrix_size/4);

% Generate the extracellular iron distribution and the corresponding field.
iron_volume = generate_iron_volume(fiber_mask_resamp,EID);
iron_field = chi_iron * generate_iron_field(iron_volume);

% Compute the total field.
field = myelin_field + iron_field;

% Calculate fiber volume fraction (FVF) and iron volume fraction (IVF).
FVF = sum(fiber_mask_resamp(:)) / (matrix_size/2)^3;
IVF = sum(iron_volume(:)) / (matrix_size/2)^3;

%%
%Calculate mGRE Signal.

% Define echo times.
TE = 0:3e-3:60e-3;
% Gyromagnetic ratio.
gamma = 2 * pi * 42.577478518e6;
% Scanner magnetic field.
B_0 = 3;
% Intraextracellular water T2.
T2_iew = 70e-3;
% Myelin water T2.
T2_myelin = 16e-3;
% Relative myelin proton density.
relative_protondensity_myelin = 0.5;

% Intra-extracellular water and myelin water weighted masks.
iew_mask = 1 - myelin_mask_resamp;
myelin_mask_resamp = relative_protondensity_myelin * myelin_mask_resamp;

% Initilize mGRE signal.
signal = zeros(1,length(TE));

% Compute the T2 decays.
T2decay_iew = exp(-TE/T2_iew);
T2decay_myelin = exp(-TE/T2_myelin);

% Computude complex factor distribution and sum it to get mGRE signal at
% TE.
for t = 1:length(TE)
complex_factor_distribution = T2decay_iew(t)*iew_mask.*exp(-1i * TE(t) * gamma * B_0 * field ) +...
                T2decay_myelin(t)*myelin_mask_resamp.*exp(-1i * TE(t) * gamma * B_0 * field );
signal(t) = sum(complex_factor_distribution(:));
end

% Normalize the signal with the first echo.
signal = signal / abs(signal(1));

% Fit for R2*.
f=fit(TE.',-log(abs(signal.')),'poly1');
R2s = f.p1;

% Save the computed parameters and the signal.
dictionary.R2s(m) = R2s;
dictionary.FVF(m) = FVF;
dictionary.g_ratio(m) = g_ratio;
dictionary.theta(m) = theta;
dictionary.IVF(m) = IVF;

dictionary.signal(m,:) = signal;

m = m+1;

end
end
end
end

% Save the global parameters to the dictionary.
dictionary.chi_iso = chi_ani;
dictionary.chi_ani = chi_iso;
dictionary.chi_iron = chi_iron;
dictionary.TE = TE;

end