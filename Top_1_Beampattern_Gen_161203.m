%
% 2016-12-03
% Sua Bae
%
% For Convex Journal
%
clc;
clear;
close all;
addpath('src');


%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transducer spec :: Convex
load('stTrans');
nDirectivityTheta = 30; % 52.5degree when width = 0.7601*lambda

% Ultrasound parameter
nSoundSpeed =  1540;
nFc = stTrans.nCenterFreq;
nWaveLength = nSoundSpeed/nFc;


asMode{1} = 'Continuous';
asMode{2} = 'Pulsed';

aAngleNum = [256,128,64,32];

asAngleDist{1} = 'Equal_alpha';
asAngleDist{2} = 'Equal_theta';

mFocalPointPos_ra = [10e-3+stTrans.nRadius, 0;...
                     10e-3+stTrans.nRadius, 15; ...
                     10e-3+stTrans.nRadius, 30; ...
                     30e-3+stTrans.nRadius, 0; ...
                     30e-3+stTrans.nRadius, 15; ...
                     30e-3+stTrans.nRadius, 30; ...   
                     30e-3+stTrans.nRadius, stTrans.nMaxTheta];     % [meter, deg]    
                 
                 
%%% Calc beampattern %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for midx = 2:-1:1 %1:2
    for anidx = 4:-1:1 %1:4
        for adidx = 1:2 %1:2
            for fidx = 7 %1:6
                
                nFolderIdx = (midx-1)*numel(asAngleDist)*numel(aAngleNum) + ...
                             anidx*numel(asAngleDist) + ... 
                             adidx;
                sFolderName = [asMode{midx}, '_Ntx', num2str(aAngleNum(anidx)), '_', asAngleDist{adidx},'\'];
                
                %%% Load plane wave angles
                load([sFolderName 'stTxAngle.mat']);
                aTxAngle    = stTxAngle.aAzi_deg;
                
                bPlotTransmitAngles = 1;
                Plot_Angles(bPlotTransmitAngles, stTxAngle);
            
                %%% TX: Aperture for each plane wave                
                mTxAperture = TxApertureForPW(stTrans, nDirectivityTheta, aTxAngle);  %%% dim: (nTransEleNum)x(nTxAngleNum)    
                mTxApod = mTxAperture;
%                 figure;imagesc(mTxAperture');

                %%% Focal point position
                aFocalPointPos_ra = mFocalPointPos_ra(fidx,:);
                [aFocalPointPos(1), aFocalPointPos(3)] = ra2xz(aFocalPointPos_ra(1),aFocalPointPos_ra(2)); % [r,a] -> [x,y,z]
                sFocalPointPos = ['r_' num2str(round((aFocalPointPos_ra(1)-stTrans.nRadius)*1e5)/1e2) '_a_' num2str(round(aFocalPointPos_ra(2)*1e2)/1e2)];
                sFocalPointPos_title = ['r=' num2str(round((aFocalPointPos_ra(1)-stTrans.nRadius)*1e5)/1e2) 'mm, a=' num2str(round(aFocalPointPos_ra(2)*1e2)/1e2) ' deg'];
                                
                %%% Find min/max available PW angle
                nTransAngle_max = stTrans.nElePitch*(stTrans.nEleNum-1)/2/stTrans.nRadius/pi*180;
                [nAvalAng_min, nAvalAng_max] = SynthAngRange(aFocalPointPos_ra(1), aFocalPointPos_ra(2), ...
                                                         -nTransAngle_max, nTransAngle_max, stTrans.nRadius, nDirectivityTheta);
                aAvalAngLogic = (nAvalAng_min <= aTxAngle)&(aTxAngle <= nAvalAng_max);
                
                %%% Available angles
                aTxAngle_trc = aTxAngle(aAvalAngLogic);
                mTxApod_trc  = mTxApod(:,aAvalAngLogic);                              
                
                
                %%% call function     
                display(['Calc ', asMode{midx}, ' Wave Beamfield Using Ntx=', num2str(aAngleNum(anidx)), '(', asAngleDist{adidx} ') at (' sFocalPointPos_title, ')']); 
                switch asMode{midx}
                    case 'Continuous'
                        [stBeamField] = CW_TX_Beampattern(stTrans, aTxAngle_trc, mTxApod_trc, aFocalPointPos, nWaveLength);
                        save([sFolderName 'stBeamField_' sFocalPointPos '.mat'], 'stBeamField');
                    case 'Pulsed'
                        [stBeamField] = PW_TX_Beampattern(stTrans, aTxAngle_trc, mTxApod_trc, aFocalPointPos, nSoundSpeed);
                        save([sFolderName 'stBeamField_' sFocalPointPos '.mat'], 'stBeamField');
                end
                
                
            end
        end
    end
end
