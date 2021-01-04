% Sua Bae
%
% 2016-12-04
%  add option, delete sDataDir, nTransEleNum, nTransElePitch, nTransRadius
%
% 2016-09-30
%  Tx beampattern/ for convex
%
% 2016-06-09
%
% mTxBeamField        : TX BeamField (consists of multiple plane waves)
% mRxBeamField        : RX BeamField for one point
%
% nWaveLength       : wave length [M] 
% aFocalPointPos    : Position of focal point [meter], dim: 2x1 (1st row: x position, 2nd row: z position)
%
% aX                : View point of x_axis [M]
% aZ                : View point of z_axis [M]
% nRxFnum             : F number (Using when bApertureGrowth is 'y')
% aPWSteerAngle     : Plane wave steering aPWSteerAngle set. [degree]
% bApertureGrowth         : 'y' - Using Aperture growth
% bApodization              : 'y' - Using RX apodization
% sApodizationWindow       : Select window (defalt : rect)
%                           'hann' - hanning window
%                           'hamm' - hamming window
% sOption          : 'TxOnly', 'RxOnly', 'TxRx'
function [mTxBeamField, mRxBeamField] = BeamFieldForOneFocalPoint_2D(   sTransType, nWaveLength, ...
                                                                        mTxSrcPos, mRxSrcPos, ...
                                                                        nTransEleWidth, ...
                                                                        aFocalPointPos, ...
                                                                        aX, aZ, aPWAngle, mTxApod_ele, aRxApod_ele, sOption)
    mTxBeamField = 0;
    mRxBeamField = 0;
    
    nXdim = length(aX);
    nZdim = length(aZ);
    
    nAngleNum = length(aPWAngle);
    
    if (strcmp(sOption,'TxOnly') || strcmp(sOption,'TxRx'))
        %% TX beamfield Calculation   
    %     mTxBeamField = 0;
        %%%% TX Beam pattern using synthetic plane wave    
        sTXType = 'PW';

        figure;
        vTxBeamField = zeros(nZdim, nXdim, nAngleNum);
        for aidx = 1:nAngleNum                
            nPWAngle = aPWAngle(aidx); % [degree]
            display(['Tx Beampattern for each PW angle ' num2str(aidx) '/' num2str(nAngleNum) ' (' num2str(nPWAngle) ' deg)']);

            %%% 1) calculate delay for a plane wave
            [aSrcDelayDist] = CalcDelayDist(sTransType, mTxSrcPos, sTXType, nPWAngle, 0, 0);

            %%%2) Apodization window dimension extension ( 1x(nTransEleNum) -> 1x(nSrcNum) )
            aTxApod_ele = mTxApod_ele(:,aidx); % (nTransEleNum)x1
            nTransEleNum = length(aTxApod_ele);
            nSrcPerEle = size(mTxSrcPos,2)/nTransEleNum; % Num. of source per element
            aTxApod_src = reshape((aTxApod_ele*ones(1,nSrcPerEle))', nTransEleNum*nSrcPerEle, 1);  % (nSrcNum)x1 (cf. EleNum x SrcPerEle = SrcNum )
            bTxDirectivity = 0;

    %         figure;
    %         for idx = 1:length(mTxSrcPos)
    %             hold on;
    %             if(aTxApod_src(idx)==0); scatter(mTxSrcPos(1,idx),mTxSrcPos(3,idx),'b');
    %             else                         scatter(mTxSrcPos(1,idx),mTxSrcPos(3,idx),'r'); end;
    %         end
    %         axis equal; set(gca,'Ydir','reverse');grid on; title('Active source (red)');  


            %%% 2) calculate a plane wave beamfield
            [vTxMag, vTxPhase]  = PressureFromEachSourceWithDelay_2D(nWaveLength, aX, aZ, mTxSrcPos, aTxApod_src, sTransType, nTransEleWidth, aSrcDelayDist, bTxDirectivity);
            vTxCmplx = vTxMag.*exp(1i*vTxPhase); % dim: (axial)x(lateral)x(source)
            mTxCmplx = sum(vTxCmplx,3); % add all beamfields from each source -> dim: (axial)x(lateral)        


            imagesc(aX*1e3,aZ*1e3,db(abs(mTxCmplx)));caxis([20 70]);colorbar; title(['Tx Beampattern ' num2str(aidx) '/' num2str(nAngleNum) ' (' num2str(nPWAngle) ' deg)']);
            pause(0.1);

            %%% 3) calculate synthetic delay for coherent compounding at focal point  = (x_f*sin(theta) + z_f*cos(theta)
            nSynDelayPhase = 2*pi/nWaveLength * (aFocalPointPos(1)*sind(nPWAngle) + aFocalPointPos(3)*cosd(nPWAngle)); % t(alpha) 
            vTxBeamField(:,:,aidx) = mTxCmplx*exp(-1i*nSynDelayPhase);   
    %         mTest = sum(vTxBeamField,3);
    %         figure;imagesc(abs(mTest));
        end

    %     save([sDataDir 'vTxBeamField'],'vTxBeamField');

        %%% 4) Synthesize (compounding) & Complex -> abs
        mTxBeamField = abs(sum(vTxBeamField,3));

    elseif (strcmp(sOption,'RxOnly') || strcmp(sOption,'TxRx'))
       %% RX beamfield Calculation

        sRXType = 'CF';    
        display(['Rx Beampattern']);

        %%% 1) Calculate delay for a plane wave
        [aRxSrcDelayDist] = CalcDelayDist(sTransType, mRxSrcPos, sRXType, 0, aFocalPointPos, 0);

        %%% 2) Apodizatio window dimension extension ( 1x(nTransEleNum) -> 1x(nSrcNum) )
        nTransEleNum = length(aRxApod_ele);
        nSrcPerEle = size(mRxSrcPos,2)/nTransEleNum; % Num. of source per element
        aRxApod_src = reshape((ones(nSrcPerEle,1)*aRxApod_ele), 1, nTransEleNum*nSrcPerEle);  % 1x(nTransEleNum) -> 1x(nSrcNum) (cf. EleNum x SrcPerEle = SrcNum )

        bRxDirectivity = 0;

        %%% 2) Calculate a RX beamfield
        [vRxMag, vRxPhase]  = PressureFromEachSourceWithDelay_2D(nWaveLength, aX, aZ, mRxSrcPos, aRxApod_src, sTransType, nTransEleWidth, aRxSrcDelayDist, bRxDirectivity);
        vRxCmplx = vRxMag.*exp(1i*vRxPhase); % dim: (axial)x(lateral)x(source)

        %%% 3) Add all beamfields from sources -> dim: (axial)x(lateral) 
        aRxCmplx = sum(vRxCmplx,3); %       

        %%% 4) Complex -> abs
        mRxBeamField = abs(aRxCmplx);
        
     end
    

end