function D_interpolated = interpolate_dictionary(D, TE_original, TE_target, degree)

if size(TE_original,2) ~= 1
    TE_original = TE_original';
end
if size(TE_target,2) ~= 1
    TE_target = TE_target';
end

A_original = [];
for i = 1:degree
    A_original = cat(2,A_original,TE_original.^(i-1));
end

A_target = [];
for i = 1:degree
    A_target = cat(2,A_target,TE_target.^(i-1));
end

D_log = log(D);

D_interpolated_log = (A_target * pinv(A_original) * D_log.').';

D_interpolated = exp(D_interpolated_log);
end