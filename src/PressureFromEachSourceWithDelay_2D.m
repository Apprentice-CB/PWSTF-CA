% 2016-06-09
% Sua Bae
%
% nWaveLength   : wavelength [meter]
% aX            : x coordinate of view points
% aZ            : z coordinate of view points
% mSrcPos       : Positions of source points dim(x,y,z,eleidx)x(SrcNum)
% aSrcAmp       : Amplitude of source points dim(1)x(SrcNum)
% aSrcDelayDist : Delay distance at each source point [meter] (dim: 1xSrcNum)
%
% vDistFromEachSrc : Distance from each source to each view point
% vMag             : Magnitude of countinuous monochromatic wave (Mag   is determined by distance from source to view point) 
% vPhase           : Phase of countinuous monochromatic wave     (Phase is determined by distance from delayed source to view point) 
%
function [vMag, vPhase] = PressureFromEachSourceWithDelay_2D(nWaveLength, aX, aZ, mSrcPos, aSrcAmp, sTransType, nTransEleWidth, aSrcDelayDist, bDirectivity)

    nXdim = length(aX);
    nZdim = length(aZ);
    nSrcNum = size(mSrcPos,2);
    
    % Increase dimension to (nZdim)x(nXdim)
    mX = repmat(aX ,[nZdim 1]); % x-coordinate of view points
    mZ = repmat(aZ',[1 nXdim]);  % z-coordinate of view points
    
    % find active source index
    aSrcIdx = (1:nSrcNum);
%     aSrcIdx = aSrcIdx(aAperture_src==1);
    
    vMag    = zeros(nZdim, nXdim, nSrcNum);
    vPhase  = zeros(nZdim, nXdim, nSrcNum);
    
    % Calculate pressure    
    for sidx = aSrcIdx

        % Distance from source
        mDistFromEachSrc = sqrt((mX - mSrcPos(1,sidx)).^2 + (mZ - mSrcPos(3,sidx)).^2);  % Distance between view positions and source points   
        switch sTransType
            case 'Linear'
                nSrc_theta = 0;
            case 'Convex'
                nSrc_theta = atan(mSrcPos(1,sidx)/mSrcPos(3,sidx)); % [rad]
        end
        switch bDirectivity
            case 0
                mBehindSrcRegion = mZ < (-tan(nSrc_theta)*(mX-mSrcPos(1,sidx)) + mSrcPos(3,sidx));
                % Mag is determined only by distance from source to point
                vMag(:,:,sidx) = aSrcAmp(sidx)./mDistFromEachSrc.*not(mBehindSrcRegion); % Amplitude of source wave
                
            case 1  
                % Directivity from source                
                mTheta = atan((mX-mSrcPos(1,sidx))./(mZ-mSrcPos(3,sidx)))-nSrc_theta;% Angle between source direction and source-to-view-point direction [rad]
                mNum = sin(pi*nTransEleWidth.*sin(mTheta)/nWaveLength);
                mDen = pi*nTransEleWidth.*sin(mTheta)/nWaveLength;
                mMul = cosd(mTheta);
                mDirectivity = mNum ./mDen .* mMul; % directivity gain
                mDirectivity(mDen==0) = 1; % NaN -> 1
                mDirectivity = max(mDirectivity,0); % minus values -> 0
                mBehindSrcRegion = mZ < (-tan(nSrc_theta)*(mX-mSrcPos(1,sidx)) + mSrcPos(3,sidx));
                mDirectivity(mBehindSrcRegion==1) = 0; % behind the element -> 0

                % Mag is determined only by distance from source to point
                vMag(:,:,sidx) = aSrcAmp(sidx)*mDirectivity./mDistFromEachSrc; % Amplitude of source wave
        end


        % Phase is determined by delay of source and distance from source to point
        nWaveNumber = 2*pi/nWaveLength;
        vPhase(:,:,sidx) = nWaveNumber*(mDistFromEachSrc + aSrcDelayDist(sidx));
    end

    
    
    
%     % Increase dimension to (nZdim)x(nXdim)x(nSrcNum)
%     vX = repmat(repmat(aX,[nZdim 1]),[1 1 nSrcNum]); % x-coordinate of view points
%     vZ = repmat(repmat(aZ',[1 nXdim]),[1 1 nSrcNum]);  % z-coordinate of view points
%     vSrcX = permute(repmat(repmat(mSrcPos(1,:),[nZdim 1]),[1 1 nXdim]),[1 3 2]); % x-coordinate of source points
%     vSrcZ = permute(repmat(repmat(mSrcPos(2,:),[nZdim 1]),[1 1 nXdim]),[1 3 2]); % z-coordinate of source points
%     vSrcDelayDist = permute(repmat(repmat(aSrcDelayDist,[nZdim 1]),[1 1 nXdim]),[1 3 2]); % z-coordinate of source points
%     
%     vDistFromEachSrc = sqrt((vX - vSrcX).^2 + (vZ - vSrcZ).^2);  % Distance between view positions and source points
%     
%     % Mag is determined only by distance from source to point
%     vMag = 1./vDistFromEachSrc; % Amplitude of source wave
%     
%     % Phase is determined by delay of source and distance from source to point
%     nWaveNumber = 2*pi/nWaveLength;
%     vPhase = nWaveNumber*(vDistFromEachSrc + vSrcDelayDist);
%     
    
    
end