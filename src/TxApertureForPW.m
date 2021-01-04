% Sua Bae 
%
% 2016-12-04
%   using stTrans
%
% 2015-09-13 
% 
% 
% Information about activated elements for each steered plane wave
%
% mTXAperture       : Logical values (1: acitve element, 0: inactive), dim = (nTransEleNum)x(AngleNum)
%
% sTransType        : 'Linear' or 'Convex'
% nTransEleNum
% nTransElePitch    : [meter]
% nTransRadius      : [meter]
% nDirectivityTheta : [degree]
% aTxAngle          : [degree] (theta, psi) when 'Matrix'
%
function [mTxAperture] = TxApertureForPW(stTrans, nDirectivityTheta, aTxAngle)

    sTransType = stTrans.sType;
    nTransEleNum = stTrans.nEleNum;
    nTransElePitch = stTrans.nElePitch;
    nTransRadius = stTrans.nRadius;
    
    
    
    switch sTransType
        
        case 'Linear'
            mTxAperture = zeros(nTransEleNum, length(aTxAngle));
            for aidx = 1:length(aTxAngle)
                mTxAperture(:,aidx) = ones(nTransEleNum,1);
            end
            
        case 'Convex'
            mTxAperture = zeros(nTransEleNum, length(aTxAngle));
            aTransEleTheta = (-(nTransEleNum-1)/2:1:(nTransEleNum-1)/2)*nTransElePitch/nTransRadius/pi*180;
            for aidx = 1:length(aTxAngle)
                
                nTxAngle = aTxAngle(aidx);  % [degree]
                
                % Half aperture size
                nTxAptrNumElHalf = round((nTransRadius*nDirectivityTheta/180*pi)/nTransElePitch); % no. of elements in 1/2 aperture.
                        
                % center element index of aperture
                nTxAptrCenterEle = round(nTransRadius*(nTxAngle-aTransEleTheta(1))/180*pi/nTransElePitch + 1);
                
                % leftmost and rightmost element index of aperture
                lft = max(nTxAptrCenterEle - nTxAptrNumElHalf, 1);
                rt  = min(nTxAptrCenterEle + nTxAptrNumElHalf, nTransEleNum);

                mTxAperture(lft:rt,aidx) = ones(rt-lft+1,1);
            end
            
        case 'Matrix'
            vTxAperture = zeros(nTransEleNum(1),nTransEleNum(2), size(aTxAngle,1));
            
            % order : (1,1), (2,1), (3,1), ..., (1,2), (2,2), (3,2), ...
            for aidx = 1:size(aTxAngle,1)
                vTxAperture(:,:,aidx) = ones(nTransEleNum(1),nTransEleNum(2),1);
            end
            mTxAperture = vTxAperture;
            
    end
    
    
    
end




