function MIMM_outputs = MIMM(dictionary,lambda_chi, QSM, Brain_Mask, iField, TE,orientation_strategy,varargin)
% Mert Sisman 7/27/2023
% 1 = Basic MIMM, 2 = Orientation Informed MIMM
if nargin < 8 &&  orientation_strategy == "orientation_informed"
    error('Theta and FA maps are required for orientation informed MIMM.');
elseif nargin == 8 &&  orientation_strategy == "orientation_informed"
    FA = varargin{1};
    theta = varargin{2};
end
switch orientation_strategy
    case "basic"
    % Basic MIMM
    lambda_theta = zeros(size(QSM));
    theta = zeros(size(QSM));
    theta_scale = 5;
    case "orientation_informed"
    % Orientation Informed MIMM. Determine the regions where DTI
    % orientation is reliable and will be used for orientation informed
    % MIMM. In other regions basic MIMM will be applied.
    FA_mask = (FA>0.25) .* Brain_Mask;
    FA_mask1 = medfilt3(FA_mask);
    FA_mask2 = imfill(FA_mask1,'holes');
    QSM_mask = QSM < 0.1;
    lambda_theta = 1e6*FA_mask2.*QSM_mask;  
    theta = theta .* Brain_Mask;
    % Round both DTI and dictionary theta values to match them exactly.
    theta_scale = 5;
    theta = round((theta)/theta_scale)*theta_scale;

end

%% Dictionary Preparation
% Interpolate the dictionary decay curves to match the data acquisition
% echo times.
dictionary.magnitude = abs(dictionary.signal); % Get dictionary magnitude curves.
TE_original =  dictionary.TE* 1e3;             % Dictionary generation TEs.
TE_target = TE * 1e3;                          % Measurement TEs.

degree = 5; % Polynomial Interpolation degree

% Interpolation
dictionary.magnitude_interpolated = interpolate_dictionary(dictionary.magnitude, TE_original, TE_target, degree);
% Normalize each dictionary element.
dictionary.magnitude_interpolated_normalized = dictionary.magnitude_interpolated./vecnorm(dictionary.magnitude_interpolated,2,2);

chi_iron = dictionary.chi_iron * 1e6; %ppm
chi_iso = dictionary.chi_iso*1e6;  %ppm

% Compute dictionary myelin volume fraction (MVF) values.
dict_MVF = dictionary.FVF.*(1-(dictionary.g_ratio).^2);
% Compute dictionary negative susceptibilities from isotropic volume
% susceptibility and MVF.
chi_neg = (chi_iso).*dict_MVF;
% Compute dictionary positive susceptibilities from isotropic volume
% susceptibility and IVF.
chi_pos = dictionary.IVF*chi_iron;
% Compute dictionary total susceptibilities that will be matched with QSM.
chi_total = chi_neg + chi_pos;
% Round dictionary theta values.
theta_dict = round((dictionary.theta)/theta_scale)*theta_scale;

%% Data Preparation
% Compute magnitude images.
magnitude = abs(iField) .* Brain_Mask;
% Normalize measurement magnitude decay curves.
magnitude = magnitude ./ vecnorm(magnitude,2,4);
magnitude(isnan(magnitude)) = 0;
magnitude(isinf(magnitude)) = 0;

% Determine the brain outer limits.
lims = bounding_box(Brain_Mask);

% Store original data size
s = size(magnitude);

% Crop the exterior regions.
Brain_Mask = Brain_Mask(lims{:});
QSM = QSM(lims{:});
magnitude = magnitude(lims{:},:).*Brain_Mask;
theta = theta(lims{:},:).*Brain_Mask;
lambda_theta = lambda_theta(lims{:},:).*Brain_Mask;

% Store cropped data size.
s0 = size(magnitude);

% Initilize original size output volumes.
MIMM_outputs.g_ratio = zeros(s(1:3));
MIMM_outputs.FVF = zeros(s(1:3));
MIMM_outputs.theta_est = zeros(s(1:3));
MIMM_outputs.R2s = zeros(s(1:3));
MIMM_outputs.chi_iron_est = zeros(s(1:3));
MIMM_outputs.chi_myelin = zeros(s(1:3));
MIMM_outputs.MVF = zeros(s(1:3));
MIMM_outputs.error = zeros(s(1:3));

% Initilize cropped size output volumes.
g_ratio0 = zeros(s0(1:3));
FVF0 = zeros(s0(1:3));
theta_est0 = zeros(s0(1:3));
R2s0 = zeros(s0(1:3));
chi_iron_est0 = zeros(s0(1:3));
chi_myelin0 = zeros(s0(1:3));
MVF0 = zeros(s0(1:3));
error0 = zeros(s0(1:3));

%% Matching Pursuit
for slice = 1:s0(3)
    % Reshape input distribution so that all the voxels in a slice can be
    % estimated simultaneously.
    mag_reshaped = reshape(magnitude(:,:,slice,:), prod(s0(1:2)), s0(4)).';
    QSM_reshaped = reshape(QSM(:,:,slice), prod(s0(1:2)), 1).';
    theta_reshaped = reshape(theta(:,:,slice), prod(s0(1:2)), 1).';
    lambda_theta_reshaped = reshape(lambda_theta(:,:,slice), prod(s0(1:2)), 1).';
    
    % Compute the error terms.
    magnitude_correlation = (abs((dictionary.magnitude_interpolated_normalized)*mag_reshaped));
    chi_error = (abs(chi_total.' - QSM_reshaped));
    theta_error = (abs(theta_dict.' - theta_reshaped));

    % Find the dictionary elements that minimize the weighted sum of the
    % error terms for each voxel.
    [error_local,ind] = min(lambda_chi*chi_error+(1-magnitude_correlation)+lambda_theta_reshaped.*theta_error);
    ind_reshaped = reshape(ind,s0(1),s0(2));
    error_local_reshaped = reshape(error_local,s0(1),s0(2));
    
    % Reshape the estimated values and insert them into the cropped output
    % volume.
    g_ratio0(:,:,slice) = dictionary.g_ratio(ind_reshaped);
    R2s0(:,:,slice) = dictionary.R2s(ind_reshaped);
    FVF0(:,:,slice) = dictionary.FVF(ind_reshaped);
    theta_est0(:,:,slice) = dictionary.theta(ind_reshaped);
    chi_iron_est0(:,:,slice) = chi_pos(ind_reshaped);
    chi_myelin0(:,:,slice) = chi_neg(ind_reshaped);
    MVF0(:,:,slice) = dict_MVF(ind_reshaped);
    error0(:,:,slice) = error_local_reshaped;
end
%% Return to Original Size by padding zeros and mask the background regions.

MIMM_outputs.g_ratio(lims{:}) = g_ratio0 .* Brain_Mask;
MIMM_outputs.FVF(lims{:}) = FVF0 .* Brain_Mask;
MIMM_outputs.theta_est(lims{:}) = theta_est0 .* Brain_Mask;
MIMM_outputs.R2s(lims{:}) = R2s0 .* Brain_Mask;
MIMM_outputs.chi_iron_est(lims{:}) = chi_iron_est0 .* Brain_Mask;
MIMM_outputs.chi_myelin(lims{:}) = chi_myelin0 .* Brain_Mask;
MIMM_outputs.MVF(lims{:}) = MVF0 .* Brain_Mask;
MIMM_outputs.error(lims{:}) = error0 .* Brain_Mask;
end


