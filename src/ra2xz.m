% Sua Bae
% 2016-09-28
%
% (radius, azimuthal angle [deg])  -> (x, z)
function [x,z] = ra2xz(r,a)
    x = r.*sind(a);
    z = r.*cosd(a);
end
