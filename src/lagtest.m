% Sua Bae
%
% 2016-12-05
%
% Pulsed wave TX beampattern using 'Field II'
%  Convex array Plane wave beamforming
%
%
function [stBeamField] = PW_TX_Beampattern(stTrans, aTxAngle_trc, mTxApod_trc_ele, aFocalPointPos_org, nSoundSpeed)
   


    bShowFieldIISetUp = 0;
    bShowPressureField = 0;
    bShowSyntheticPattern = 0;
        
    %%
    aFocalPointPos(1) = aFocalPointPos_org(1);
    aFocalPointPos(2) = aFocalPointPos_org(2);
    aFocalPointPos(3) = aFocalPointPos_org(3) - stTrans.nRadius;
    
    %% View points
    aX = (-40: 1/5 : 40)*1e-3 + aFocalPointPos(1);
    aZ = (-15: 1/5 : 15)*1e-3 + aFocalPointPos(3);    
    [vX, vY, vZ]= ndgrid(aX,0,aZ);
    vPoints(:,1) = vX(:);
    vPoints(:,2) = vY(:);
    vPoints(:,3) = vZ(:);
    
    %% Field II setup
    addpath('C:\Field_II');
    
    %%% Transducer spec 
    stTrans.sType       = 'Convex';
    stTrans.CenterFreq  = 3.125e6;
    stTrans.nRadius     = 41.219e-3;
    stTrans.nEleNum_x   = stTrans.nEleNum;
    stTrans.nEleNum_y   = 1;
    stTrans.nEleNum     = stTrans.nEleNum_x * stTrans.nEleNum_y;
    stTrans.nPitch_x    = stTrans.nElePitch;
    stTrans.nPitch_y    = 0;
    stTrans.nWidth_x    = stTrans.nEleWidth;
    stTrans.nWidth_y    = 6e-3;
    
    %%% Ultrasound parameter
    nAttenFactor 	= 0.5/(1e-2*1e6);   % Attenuation factor[dB/(m*Hz)] % 0.5 dB/(cm*MHz)
    nFc             = stTrans.CenterFreq; % [Hz]
    nFs             = stTrans.CenterFreq * 4; % [Hz]
    nUpSampFactor   = 8;
    nSimulFrequency	= nFs*nUpSampFactor; % [Hz]
    nWaveLength 	= nSoundSpeed/nFc;    
    dt              = 1/nSimulFrequency; % [sec]
    nTc             = 1/nFc;        % period [sec]
    nTxPulseCycle   = 1;
    
    %%% Field
    field_init(0);
    set_field('c',nSoundSpeed);             % Set speed of sound
    set_field('Freq_att',nAttenFactor);     % frequency dependency attenuation in [dB/(m*Hz)] around the center frequency
    set_field('att', nFc*nAttenFactor);     % Freqency independent attenuation[dB/m] 
    set_field('att_f0',nFc);                % Attenuation center frequency[Hz]
    set_field('use_att',1);                 % attenuation on/off
    set_field('fs',nSimulFrequency);        % set the sampling frequency
    set_field('use_rectangles',1);          % use ractangles for apertures

    %%% Transducer
    nFarFieldDepth_x = 0.01e-3; % depth over nFarFieldDepth_x is assumed to be far field on x-z plane 
    nFarFieldDepth_y = 0.1e-3; % depth over nFarFieldDepth_x is assumed to be far field on y-z plane 

    nMathEleSize_x = sqrt(nFarFieldDepth_x * 4 * nWaveLength); % H < sqrt(4*lambda*z)
    nMathEleSize_y = sqrt(nFarFieldDepth_y * 4 * nWaveLength);

    nSubDivNum_x = ceil(stTrans.nWidth_x / nMathEleSize_x);
    nSubDivNum_y = ceil(stTrans.nWidth_y/ nMathEleSize_y);

    % Has no elevational focusing lens 
    pTxTransducer = xdc_convex_array (stTrans.nEleNum_x, ...
                                      stTrans.nWidth_x, ...
                                      stTrans.nWidth_y, ...
                                      stTrans.nPitch_x-stTrans.nWidth_x, ...
                                      stTrans.nRadius, ...
                                      nSubDivNum_x, ...
                                      nSubDivNum_y, ...
                                      aFocalPointPos);
    if(bShowFieldIISetUp)
        figure; 
        show_xdc(pTxTransducer); title('Tx Aperture');
    end
    
    %  ??
    xdc_baffle(pTxTransducer,0);      

    %%% Impulse Response (-6dB BW is 75% when nFc = 3.125 MHz)
    at2 = 0:dt:2.8*nTc;
    aTransImpulsResp = sin(2*pi*nFc*at2); 
    aTransImpulsResp = aTransImpulsResp.*(hanning(max(size(aTransImpulsResp)))'); % Transducer's impulse response
    xdc_impulse(pTxTransducer,aTransImpulsResp); % Tx impulse response 
        
    %%% Excitation Pulse
    at1 = 0:dt:nTxPulseCycle*nTc;
    nPulseAmp = 1;
    aPulseSeq  = nPulseAmp*sin(2*pi*nFc*at1);  % excitation pulse
    xdc_excitation(pTxTransducer,aPulseSeq);  % convol Tx aperture and pulse

    if(bShowFieldIISetUp)
        % Excitation pulse, Impulse response
        figure;
        subplot(3,1,1)
        plot(at1*1e6,aPulseSeq), title('Exitation Pulse')
        subplot(3,1,2)
        plot(at2*1e6,aTransImpulsResp), title('Transducer`s Impulse response')
        subplot(3,1,3)
        aConvSeq = conv(conv(aPulseSeq,aTransImpulsResp),aTransImpulsResp);
        plot(0:dt*1e6:dt*(length(aConvSeq)-1)*1e6,aConvSeq), title('Expected receive pulse')
    end
%     nLag = round(((numel(aPulseSeq)+numel(aTransImpulsResp)-1)+numel(aTransImpulsResp)-1)/2); % for received signal
    nLag = round((numel(aPulseSeq)+numel(aTransImpulsResp)-1)/2); % for pressure field
        
    %%% Take back source information for delay calculation
    nSrcNum = nSubDivNum_x*nSubDivNum_y;
    mSrcPos = GetSrcPos(pTxTransducer, stTrans.nEleNum_x, stTrans.nEleNum_y, nSubDivNum_x, nSubDivNum_y, stTrans.sType, stTrans.nRadius);
    
    
    %% Tx beampattern for every plane wave
    
    nSample  = round((aFocalPointPos_org(3)+50e-3)/nSoundSpeed*nFs);
    mP_syn = zeros(nSample,size(vPoints,1));
    for pwidx = 1:numel(aTxAngle_trc)
        
        %%% Angle
        nTxAngle = aTxAngle_trc(pwidx);
        display(['Tx Beampattern for each PW angle ' num2str(pwidx) '/' num2str(numel(aTxAngle_trc)) ' (' num2str(nTxAngle) ' deg)']);
        
        %%% Apodization
        aTxApod = mTxApod_trc_ele(:,pwidx);
        xdc_apodization(pTxTransducer, 0, aTxApod');
        
        %%% Delay
        aTxDelayDist_src = CalcDelayDist(stTrans.sType, mSrcPos, 'PW', nTxAngle, 0, 0);
        aTxDelayDist_ele = aTxDelayDist_src(1:nSrcNum:end);
        aTxDelay =  aTxDelayDist_ele/nSoundSpeed; % [sec] <- [meter]
        xdc_focus_times(pTxTransducer,0,aTxDelay); % Tx Delay setting
        
        
        %%% Calculate the received signals from a collection of scatterers for all elements
        [mP_org, t] = calc_hp(pTxTransducer, vPoints); 

        %%% Zero-padding (so that the first sample is captured at t = 0)
        nPadLen = round(t/dt);
        if (nPadLen < 0)
            mP_pad = mP_org(abs(nPadLen):end,:);
        else
            mP_pad = padarray(mP_org,[nPadLen 0],'pre');
        end
        clear mP_org

        %%% Truncation or Zero-Padding (so that # of samples = RFInfo.nSample)
        nUpSample = nSample*nUpSampFactor;
        if (size(mP_pad,1) > nUpSample)
            mP_trc = mP_pad(1:nUpSample,:);
        elseif(size(mP_pad,1) < nUpSample)
            mP_trc = padarray(mP_pad,[nUpSample-size(mP_pad,1) 0],'post');
        end
        clear mP_pad 

        %%% Calculate synthetic delay for coherent compounding
        nSynDist = aFocalPointPos_org(1)*sind(nTxAngle) + aFocalPointPos_org(3)*cosd(nTxAngle);
        nSynDelayOffset_px = round( (nSynDist-aFocalPointPos_org(3))/nSoundSpeed*nSimulFrequency )  +  nLag; % [pixel]

        if (nSynDelayOffset_px > 0)
            mP_trc = padarray(mP_trc(nSynDelayOffset_px+1:end,:),[nSynDelayOffset_px 0],'post');            
        else
            mP_trc = padarray(mP_trc(1:end-abs(nSynDelayOffset_px),:),[abs(nSynDelayOffset_px) 0],'pre');
        end

        %%% Decimate to lower smapling frequency (nFs)
        mP_deci = mP_trc(1:nUpSampFactor:end,:);
        clear mP_trc
        
        if(bShowPressureField) 
            vP_deci = permute(reshape(mP_deci,nSample,numel(aX),numel(aZ)),[3,2,1]);
            figure;
            for idx = 1:nSample
                imagesc(aX*1e3, aZ*1e3, squeeze(vP_deci(:,:,idx)));title(num2str(idx));
                axis equal; axis tight; xlabel('x [mm]'); ylabel('z [mm]');
                pause(0.001);
            end
        end
        
        mP_syn = mP_syn + mP_deci; 
        clear mP_deci
        
        if(bShowSyntheticPattern)
            vP = permute(reshape(mP_syn,nSample,numel(aX),numel(aZ)),[3,2,1]);
            figure;
            imagesc(aX*1e3, aZ*1e3, squeeze(vP(:,:,round(aFocalPointPos_org(3)/nSoundSpeed*nFs))));title(['t=' num2str(round(aFocalPointPos_org(3)/nSoundSpeed*1e3*100)/100) ' ms']);;
            axis equal; axis tight; xlabel('x [mm]'); ylabel('z [mm]');
            pause(0.001);            
%             for idx = 300:nSample
%                 imagesc(aX*1e3, aZ*1e3, squeeze(vP(:,:,idx)));title(num2str(idx));
%                 axis equal; axis tight; xlabel('x [mm]'); ylabel('z [mm]');
%                 pause(0.001);
%             end
        end
        
        
    end
    
    mIntensity = mP_syn.*mP_syn;
    clear mtvP
    
    vIntensity_t = permute(reshape(mIntensity, nSample, numel(aX), numel(aZ)),[2,3,1]);
    
    %% Export
    stBeamField.aFocalPointPos      = aFocalPointPos;
    stBeamField.aX                  = aX;
    stBeamField.aZ                  = aZ;
    stBeamField.vIntensity_t        = vIntensity_t;
    stBeamField.aTxAngle_trc        = aTxAngle_trc;
    stBeamField.mTxApod_trc_ele     = mTxApod_trc_ele;

    
end
    



