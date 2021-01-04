% Sua Bae
%
% 2016-12-04
%   correction 
% 2016-10-26
%   calculate minimum and maximum angle available for PWSF
%
% aP_r: radius of focal point [m]
% aP_a: angle of focal point [deg]
% nTransAngle_min: minimum angle of txdcr elements [deg]
% nTransAngle_max: maximum angle of txdcr elements [deg]
% nTransRadius: radius of txdcr [m]
% nTransAccAngle: acceptance angle of txdcr [deg]
%
function [aTheta_min, aTheta_max] = SynthAngRange(aP_r, aP_a, nTransAngle_min, nTransAngle_max, nTransRadius, nTransAccAngle)

    aP_x = aP_r.*sind(aP_a);
    aP_z = aP_r.*cosd(aP_a);
    
    nTransLft_x = nTransRadius*sind(nTransAngle_min);
    nTransLft_z = nTransRadius*cosd(nTransAngle_min);
    nTransRgt_x = nTransRadius*sind(nTransAngle_max);
    nTransRgt_z = nTransRadius*cosd(nTransAngle_max);

    % 1. Maximum range of synthetic plane wave angle (neglecting directivity angle)
    aAng1_max = atand((aP_x-nTransLft_x)./(aP_z-nTransLft_z));
    aAng1_min = atand((aP_x-nTransRgt_x)./(aP_z-nTransRgt_z));
    
    % 2. Maximum range of synthetic plane wave angle (considering directivity angle)
    aAng2_max_tmp = aP_a + asind((nTransRadius*sind(nTransAccAngle))./aP_r);
    aAng2_min_tmp = aP_a - asind((nTransRadius*sind(nTransAccAngle))./aP_r);
    aAng2_max = min(aAng2_max_tmp, aAng1_max);
    aAng2_min = max(aAng2_min_tmp, aAng1_min);
    
    
    aTheta_max = aAng2_max;
    aTheta_min = aAng2_min;
    
    
end