function [iron_volume] = generate_iron_volume(fiber_volume,desired_iron_density)
s = size(fiber_volume);
iron_volume = (rand(s)<desired_iron_density) .* (1-fiber_volume);
end