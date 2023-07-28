function [iron_field] = generate_iron_field(iron_volume)
s = size(iron_volume);

D = dipole_kernel(s*2, [1,1,1], [0,0,1]);

susceptibility_3D = padarray(iron_volume,s/2,0,'both');

iron_field = real(ifftn(fftn(susceptibility_3D).*D));

iron_field = iron_field(s(1)/2+1:end-s(1)/2,s(2)/2+1:end-s(2)/2,s(3)/2+1:end-s(3)/2);
end