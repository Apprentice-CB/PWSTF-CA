% Sua Bae
%
% 2016-12-04
%
% Continuous wave TX beampattern using 'BeamFieldForOneFocalPoint_2D.m'
%  Convex array Plane wave beamforming
%
function [stBeamField] = CW_TX_Beampattern(stTrans, aTxAngle_trc, mTxApod_trc_ele, aFocalPointPos, nWaveLength)


    sTransType = stTrans.sType;
    nTransEleNum = stTrans.nEleNum;
    nTransElePitch = stTrans.nElePitch;
    nTransEleWidth = stTrans.nEleWidth;
    nTransRadius = stTrans.nRadius;
    
    %% TX: Generate Sources

    nSrcNumPerTransEle = 5; % Number of source points per element
    
    nTxAptrCenterPos = 0;
    nTxAptrCenterAngle = 0;  

    % source poisition and element index, 
    mTxSrcPos = GenerateSource(sTransType, nTransElePitch, nTransEleWidth, nTransRadius, ...
                                nTransEleNum, nSrcNumPerTransEle, nTxAptrCenterPos, nTxAptrCenterAngle); % dim = (xpos, zpos, eleidx)x(nRxAptrEleNum*nSrcNumPerTransEle)
    hold on;
    scatter(aFocalPointPos(1), aFocalPointPos(3));

    %% View points
    aX = (-40: 1/10 : 40)*1e-3 + aFocalPointPos(1);
    aZ = (-15: 1/10 : 15)*1e-3 + aFocalPointPos(3);
%     aZ = (-5: 1/10 : 5)*1e-3 + aFocalPointPos(3);
%     aZ = aFocalPointPos(3);

    %% Simulation

    [mTxBeamField, tmp] = BeamFieldForOneFocalPoint_2D(sTransType, nWaveLength, ...
                                                       mTxSrcPos, 0, ... % (x,y,z,eidx)
                                                       nTransEleWidth, ...
                                                       aFocalPointPos, ... % (x,y,z)
                                                       aX, aZ, aTxAngle_trc, mTxApod_trc_ele,  0, 'TxOnly');


    %% Export
    stBeamField.aFocalPointPos      = aFocalPointPos;
    stBeamField.aX                  = aX;
    stBeamField.aZ                  = aZ;
    stBeamField.mTxBeamField        = mTxBeamField;
    stBeamField.aTxAngle_trc        = aTxAngle_trc;
    stBeamField.mTxApod_trc_ele     = mTxApod_trc_ele;

    
end
    



