function [field] = compute_field(matrix_size,r_o,g_ratio, theta, chi_iso, chi_ani)
%% Compute 2D field dirstribution given hollow circular distribution.
% Derived in Wharton, S. and Bowtell, R., 2012. Fiber orientation-dependent white matter contrast in gradient echo MRI.
% Proceedings of the National Academy of Sciences, 109(45), pp.18559-18564.
% Output: field = total field perturbation.
% Inputs: matrix_size = matrix size
%         r_o = outer radius including myelin.
%         g_ratio = r_i / r_o
%         theta = fiber orientation.
%         chi_iso = isotropic susceptibility of myelin.
%         chi_ani = anisotropic susceptibility of myelin.

%%
r_i = ceil(r_o*g_ratio); % inner radius excluding myelin.
% Generate cartesion grid.
lim = (matrix_size - 1) / 2;
x = -lim:lim;
y = -lim:lim;
% Convert it to polar grid.
[phi,r] = polar_mesh(x,y);

% Create masks that will respresent inner (intracellular), annular (myelin)
% and outer (extracellular) regions.
inner_region = double(abs(r) < r_i);
annular_region = double((abs(r) >= r_i) .* (abs(r) < r_o));
outer_region = double(abs(r) >= r_o);

%%
%Field due to isotropic susceptibility in 3 regions.
delta_f_outer_iso = (chi_iso .* sind(theta)^2 .* cosd(2*phi)) ./ 2 .* ((r_o^2 - r_i^2)./(r).^2);
delta_f_annular_iso = chi_iso ./ 2 .* (cosd(theta)^2 - 1/3 - sind(theta)^2 .* cosd(2*phi).* (r_i^2./r.^2));
delta_f_inner_iso = 0;

delta_f_iso = delta_f_inner_iso .* inner_region + delta_f_annular_iso .* annular_region + delta_f_outer_iso .* outer_region;


%%
%Field due to anisotropic susceptibility in 3 regions.

delta_f_outer_ani = (chi_ani .* (sind(theta)).^2 .* cosd(2*phi))./8 .* ((r_o^2 - r_i^2)./r.^2);
delta_f_annular_ani = chi_ani .* (sind(theta).^2.*(-5/12 - cosd(2*phi)/8 .* (1+r_i^2./r.^2) + 3/4 .* log(r_o./r))-cosd(theta)^2/6);
delta_f_inner_ani = (3 .* chi_ani  .* sind(theta).^2) / 4 .* log(r_o./r_i);


delta_f_ani = delta_f_inner_ani .* inner_region + delta_f_annular_ani .* annular_region + delta_f_outer_ani .* outer_region;
delta_f_ani(isnan(delta_f_ani)) = 0;


%%
%Total field.

field = delta_f_ani + delta_f_iso;
end
