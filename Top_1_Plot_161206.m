% Sua Bae
%
% 2016-12-06
%  Comparing the beamfield: Equal alpha vs. equal theta
%    1) CW simul
%    2) PW simul

clear;
close all;

addpath('src');

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)   
%% View option

bCWsimul = 0;
bPWsimul = 1;


%% Load Transducer 
load('stTrans');
nTransAccAngle = 30;

%% Data set params
asMode{1} = 'Continuous';
asMode{2} = 'Pulsed';
aAngleNum = [256,128,64,32];
asAngleDist{1} = 'Equal_alpha';
asAngleDist{2} = 'Equal_theta';
mFocusPos_ra = [10e-3+stTrans.nRadius, 0;...
                10e-3+stTrans.nRadius, 15; ...
                10e-3+stTrans.nRadius, 30; ...
                30e-3+stTrans.nRadius, 0; ...
                30e-3+stTrans.nRadius, 15; ...
                30e-3+stTrans.nRadius, 30];     % [meter, deg]     
%% 1. CW simul

if(bCWsimul)
    midx = 1; % continuous wave

    %   1. 1) Varying Ntx, focus

    for anidx = 4   
        for fidx = 4:6%1:6

            %%% Focus
            aFocusPos_ra = mFocusPos_ra(fidx,:);
            [aFocusPos(1), aFocusPos(3)] = ra2xz(aFocusPos_ra(1),aFocusPos_ra(2)); % [r,a] -> [x,y,z]
            sFocusPos = ['r_' num2str(round((aFocusPos_ra(1)-stTrans.nRadius)*1e5)/1e2) '_a_' num2str(round(aFocusPos_ra(2)*1e2)/1e2)];
            sFocusPos_title = ['r=' num2str(round((aFocusPos_ra(1)-stTrans.nRadius)*1e5)/1e2) 'mm, a=' num2str(round(aFocusPos_ra(2)*1e2)/1e2) ' deg'];


            %%% Load Data
            adidx = 1;
                %%% Folder
                sFolderName = [asMode{midx}, '_Ntx', num2str(aAngleNum(anidx)), '_', asAngleDist{adidx},'\'];
                %%% Load data
                load([sFolderName 'stBeamField_' sFocusPos '.mat']);
                stBF_ea         = stBeamField;
                load([sFolderName 'stTxAngle.mat']);
                aTxAngle_ea   = stTxAngle.aAzi_deg;
                
            adidx = 2;
                %%% Folder
                sFolderName = [asMode{midx}, '_Ntx', num2str(aAngleNum(anidx)), '_', asAngleDist{adidx},'\'];
                %%% Load data
                load([sFolderName 'stBeamField_' sFocusPos '.mat']);
                stBF_et         = stBeamField;
                load([sFolderName 'stTxAngle.mat']);
                aTxAngle_et   = stTxAngle.aAzi_deg;

            %%% Plot
            figure('Position',[500, 100, 600, 800]);   
            aX              = stBF_ea.aX;
            aZ              = stBF_ea.aZ;
            mIntensity_ea   = stBF_ea.mTxBeamField;
            mIntensity_et   = stBF_et.mTxBeamField;
            [mZ, mX] = ndgrid(aZ', aX');
            aIntensity_ea = interpn(mZ, mX, mIntensity_ea, aFocusPos(3), aX);
            aIntensity_et = interpn(mZ, mX, mIntensity_et, aFocusPos(3), aX);
            subplot(3,1,1); % xz plane view           
                imagesc(aX*1e3, (aZ-stTrans.nRadius)*1e3, db(mIntensity_ea)); 
                axis tight; axis equal; xlabel('x [mm]'); zlabel('z [mm]'); caxis([db(max(aIntensity_ea))-60, db(max(aIntensity_ea))]);
                title('Equal {\Delta}{\alpha}'); colorbar;
            subplot(3,1,2); % xz plane view 
                imagesc(aX*1e3, (aZ-stTrans.nRadius)*1e3, db(mIntensity_et));   
                axis tight; axis equal; xlabel('x [mm]'); zlabel('z [mm]'); caxis([db(max(aIntensity_et))-60, db(max(aIntensity_ea))]);
                title('Equal {\Delta}{\theta}'); colorbar;
            subplot(3,1,3); % x-axis view            
                plot(aX*1e3, db(aIntensity_ea),'LineWidth',1.5); 
                hold on;
                plot(aX*1e3, db(aIntensity_et),'LineWidth',1.5);    
                legend('Equal {\Delta}{\alpha}','Equal {\Delta}{\theta}','Location','NorthEast');
                title('Beamfield on x-axis [dB]'); 
                ylim([db(max(aIntensity_ea))-50, db(max(aIntensity_ea))+10]);
                xlim([aX(1)*1e3, aX(end)*1e3]); xlabel('x-axis [mm]');
                grid on; grid minor;

%             %%% Plot angle set
%             figure('Position',[100, 600, 400, 200]); 
%             aTxAngle_trc_ea   = stBF_ea.aTxAngle_trc;
%             aTxAngle_trc_et   = stBF_et.aTxAngle_trc;
%             subplot(2,1,1);      
%                 stem(aTxAngle_trc_ea, ones(numel(aTxAngle_trc_ea),1),'filled','LineWidth',1.5, 'MarkerSize', 3);
%                 grid minor; set(gca, 'XTick', -40:20:40); set(gca, 'YTick', []); xlim([-52,52]);
%             subplot(2,1,2);   
%                 stem(0,0,'MarkerSize', 0.1); hold on;
%                 stem(aTxAngle_trc_et, ones(numel(aTxAngle_trc_et),1),'filled','LineWidth',1.5, 'MarkerSize', 3);            
%                 grid minor; set(gca, 'XTick', -40:20:40); set(gca, 'YTick', []); xlim([-52,52]);
% 
%             %%% Plot SIN OF angle set
%             figure('Position',[100, 200, 400, 200]); 
%             aTxAngle_trc_ea   = stBF_ea.aTxAngle_trc;
%             aTxAngle_trc_et   = stBF_et.aTxAngle_trc;
%             subplot(2,1,1);      
%                 stem(sind(aTxAngle_trc_ea), ones(numel(aTxAngle_trc_ea),1),'filled','LineWidth',1.5, 'MarkerSize', 3);
%                 grid minor; set(gca, 'XTick', -0.6:0.2:0.6); set(gca, 'YTick', []); xlim([sind(-52),sind(52)]);
%             subplot(2,1,2);   
%                 stem(0,0,'MarkerSize', 0.1); hold on;
%                 stem(sind(aTxAngle_trc_et), ones(numel(aTxAngle_trc_et),1),'filled','LineWidth',1.5, 'MarkerSize', 3);            
%                 grid minor; set(gca, 'XTick', -0.6:0.2:0.6); set(gca, 'YTick', []); xlim([sind(-52),sind(52)]);

        end
    end
end

%% 2. PW simul
if(bPWsimul)

    midx = 2; % pulsed wave

    %   1. 1) Varying Ntx, focus

    for anidx = 4   
        for fidx = 4%1:6

            %%% Focus
            aFocusPos_ra = mFocusPos_ra(fidx,:);
            [aFocusPos(1), aFocusPos(3)] = ra2xz(aFocusPos_ra(1),aFocusPos_ra(2)); % [r,a] -> [x,y,z]
            sFocusPos = ['r_' num2str(round((aFocusPos_ra(1)-stTrans.nRadius)*1e5)/1e2) '_a_' num2str(round(aFocusPos_ra(2)*1e2)/1e2)];
            sFocusPos_title = ['r=' num2str(round((aFocusPos_ra(1)-stTrans.nRadius)*1e5)/1e2) 'mm, a=' num2str(round(aFocusPos_ra(2)*1e2)/1e2) ' deg'];


            %%% Load Data
            adidx = 1;
                %%% Folder
                sFolderName = [asMode{midx}, '_Ntx', num2str(aAngleNum(anidx)), '_', asAngleDist{adidx},'\'];
                %%% Load data
                load([sFolderName 'stBeamField_' sFocusPos '.mat']);
                stBF_ea         = stBeamField;                
                load([sFolderName 'stTxAngle.mat']);
                aTxAngle_ea   = stTxAngle.aAzi_deg;
            adidx = 2;
                %%% Folder
                sFolderName = [asMode{midx}, '_Ntx', num2str(aAngleNum(anidx)), '_', asAngleDist{adidx},'\'];
                %%% Load data
                load([sFolderName 'stBeamField_' sFocusPos '.mat']);
                stBF_et         = stBeamField;       
                load([sFolderName 'stTxAngle.mat']);
                aTxAngle_et   = stTxAngle.aAzi_deg;

            %%% Plot x-axis view
            figure('Position',[500, 100, 600, 800]);   
            aX              = stBF_ea.aX;
            aZ              = stBF_ea.aZ + stTrans.nRadius;
            mIntensity_ea   = stBF_ea.vIntensity_t(:,:,stBF_ea.nTime_at_Focus_px) * 1e26;
            mIntensity_et   = stBF_et.vIntensity_t(:,:,stBF_et.nTime_at_Focus_px) * 1e26;
            [mZ, mX] = ndgrid(aZ', aX');
            aIntensity_ea = interpn(mZ, mX, mIntensity_ea, aFocusPos(3), aX);
            aIntensity_et = interpn(mZ, mX, mIntensity_et, aFocusPos(3), aX);            
            aLine_x = [aX(1), aX(end)]*1e3;
            aLine_z = [aFocusPos(3)*1e3, aFocusPos(3)*1e3] - stTrans.nRadius*1e3;
            subplot(3,1,1);            
                imagesc(aX*1e3, (aZ-stTrans.nRadius)*1e3, db(mIntensity_ea)); 
                hold on;line(aLine_x,aLine_z,'Color','w','LineStyle','--');
                axis tight; axis equal; xlabel('x [mm]'); zlabel('z [mm]'); caxis([db(max(aIntensity_ea))-100, db(max(aIntensity_ea))]);
                title(['Equal {\Delta}{\alpha}, (Ntx=' num2str(aAngleNum(anidx)),')']); colorbar;
            subplot(3,1,2);   
                imagesc(aX*1e3, (aZ-stTrans.nRadius)*1e3, db(mIntensity_et));   
                hold on;line(aLine_x,aLine_z,'Color','w','LineStyle','--');
                axis tight; axis equal; xlabel('x [mm]'); zlabel('z [mm]'); caxis([db(max(aIntensity_et))-100, db(max(aIntensity_et))]);
                title(['Equal {\Delta}{\theta}, (Ntx=' num2str(aAngleNum(anidx)),')']); colorbar;
            subplot(3,1,3);
                plot(aX*1e3, db(aIntensity_ea),'LineWidth',1.5); 
                hold on;
                plot(aX*1e3, db(aIntensity_et),'LineWidth',1.5);    
                legend('Equal {\Delta}{\alpha}','Equal {\Delta}{\theta}','Location','NorthEast');
                title('Beamfield on x-axis [dB]');                 
                ylim([db(max(aIntensity_ea))-100, db(max(aIntensity_ea))+10]);
                xlim([aX(1)*1e3, aX(end)*1e3]); xlabel('x-axis [mm]');
                grid on; grid minor;
                
            %%% Plot SIN OF angle set with shade!
            figure('Position',[100, 200, 400, 200]);
            [nMinTheta, nMaxTheta] = SynthAngRange(aFocusPos_ra(1), aFocusPos_ra(2), min(aTxAngle_ea), max(aTxAngle_ea), stTrans.nRadius, nTransAccAngle);
            subplot(2,1,1);      
                rectangle('Position',[sind(nMinTheta), 0, sind(nMaxTheta)-sind(nMinTheta), 1], 'FaceColor', [0.8,0.8,0.8], 'EdgeColor', [0.8,0.8,0.8]);
                hold on; stem(sind(aTxAngle_ea), ones(numel(aTxAngle_ea),1),'filled','LineWidth',1.5, 'MarkerSize', 3);
                grid minor; set(gca, 'XTick', -0.6:0.2:0.6); set(gca, 'YTick', []); xlim([sind(-52),sind(52)]);
            subplot(2,1,2); 
                rectangle('Position',[sind(nMinTheta), 0, sind(nMaxTheta)-sind(nMinTheta), 1], 'FaceColor', [0.8,0.8,0.8], 'EdgeColor', [0.8,0.8,0.8]);
                hold on; stem(0,0,'MarkerSize', 0.1); 
                hold on; stem(sind(aTxAngle_et), ones(numel(aTxAngle_et),1),'filled','LineWidth',1.5, 'MarkerSize', 3);            
                grid minor; set(gca, 'XTick', -0.6:0.2:0.6); set(gca, 'YTick', []); xlim([sind(-52),sind(52)]);
                
%             %%% Plot angle set
%             figure('Position',[100, 600, 400, 200]); 
%             aTxAngle_trc_ea   = stBF_ea.aTxAngle_trc;
%             aTxAngle_trc_et   = stBF_et.aTxAngle_trc;
%             subplot(2,1,1);      
%                 stem(aTxAngle_trc_ea, ones(numel(aTxAngle_trc_ea),1),'filled','LineWidth',1.5, 'MarkerSize', 3);
%                 grid minor; set(gca, 'XTick', -40:20:40); set(gca, 'YTick', []); xlim([-52,52]);
%             subplot(2,1,2);   
%                 stem(0,0,'MarkerSize', 0.1); hold on;
%                 stem(aTxAngle_trc_et, ones(numel(aTxAngle_trc_et),1),'filled','LineWidth',1.5, 'MarkerSize', 3);            
%                 grid minor; set(gca, 'XTick', -40:20:40); set(gca, 'YTick', []); xlim([-52,52]);
% 
%             %%% Plot SIN OF angle set
%             figure('Position',[100, 200, 400, 200]); 
%             aTxAngle_trc_ea   = stBF_ea.aTxAngle_trc;
%             aTxAngle_trc_et   = stBF_et.aTxAngle_trc;
%             subplot(2,1,1);      
%                 stem(sind(aTxAngle_trc_ea), ones(numel(aTxAngle_trc_ea),1),'filled','LineWidth',1.5, 'MarkerSize', 3);
%                 grid minor; set(gca, 'XTick', -0.6:0.2:0.6); set(gca, 'YTick', []); xlim([sind(-52),sind(52)]);
%             subplot(2,1,2);   
%                 stem(0,0,'MarkerSize', 0.1); hold on;
%                 stem(sind(aTxAngle_trc_et), ones(numel(aTxAngle_trc_et),1),'filled','LineWidth',1.5, 'MarkerSize', 3);            
%                 grid minor; set(gca, 'XTick', -0.6:0.2:0.6); set(gca, 'YTick', []); xlim([sind(-52),sind(52)]);

        end
    end
end

