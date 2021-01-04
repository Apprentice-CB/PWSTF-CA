% Sua Bae
%
% 2016-11-04
%  Summation of exponentionals..
%
% mX: x position of observation field
% mZ: z position of observation field
% aFocalPointPos: (x,y,z) position of focus
function mTxBeamField_sumexp = NumericBeamPattern_SumOfExp(mX, mZ, aFocalPointPos, aTxAngle_trc, nWaveLength)


    mComplex_env = zeros(size(mX));
    
    for idx = 1:numel(aTxAngle_trc)
        nDelayDist = aFocalPointPos(1)*sind(aTxAngle_trc(idx)) + aFocalPointPos(3)*cosd(aTxAngle_trc(idx));
        mComplex_env = mComplex_env + exp(1i*2*pi/nWaveLength*(mX*sind(aTxAngle_trc(idx))+mZ*cosd(aTxAngle_trc(idx)) - nDelayDist));
    end
    
    mTxBeamField_sumexp = abs(mComplex_env);
    
end