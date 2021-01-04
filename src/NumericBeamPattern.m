% Sua Bae
%
% 2016-11-03
%  Add 'bUsingMontaldoApprox'
%
% 2016-10-28
%
% Numerical model of beampattern
%
% aX                    : observational axis
% nAxRotAngle           : When the observational axes are rotated, PW angles are also rotated
% aFocalPointPos        : Focal point position
% aTxAngle_trc          : PW angles used for PWSTF
% nWavelength           : wavelength
% bUsingMontaldoApporx  : Delta alpha in the sinc function is replaced with delta theta instead of average delta alpha
function aTxBeamField_numeric = NumericBeamPattern(aX, nAxRotAngle, aFocalPointPos, aTxAngle_trc, nWaveLength, bUsingMontaldoApprox)

    aX_sh = aX - aFocalPointPos(1);    
    nTxAngleNum_trc = numel(aTxAngle_trc);
    
    if (bUsingMontaldoApprox)
        dTheta = aTxAngle_trc(2)-aTxAngle_trc(1);
        dTxAlpha = sind(dTheta);
    else
        dTxAlpha = mean( sind(aTxAngle_trc(2:end)-nAxRotAngle) - sind(aTxAngle_trc(1:end-1)-nAxRotAngle) );
    end

    aNum = sin(pi*aX_sh*dTxAlpha*nTxAngleNum_trc/nWaveLength);
    aDen = sin(pi*aX_sh*dTxAlpha/nWaveLength);    
    aTxBeamField_numeric = aNum./aDen;
    aTxBeamField_numeric(aDen==0) = nTxAngleNum_trc;
    
end