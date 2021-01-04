%
% 2016-12-03
% Sua Bae
%
% For Convex Journal
%
clc;
clear;
close all;


%% Parameters

asMode{1} = 'Continuous';
asMode{2} = 'Pulsed';

aAngleNum = [256,128,64,32];

asAngleDist{1} = 'Equal_alpha';
asAngleDist{2} = 'Equal_theta';

mFocus = [10e-3, 0; ...
          10e-3, 15; ...
          10e-3, 30; ...
          30e-3, 0; ...
          30e-3, 15; ...
          30e-3, 30];          

%% Angle Gen 

for midx = 1:2
    for anidx = 1:4
        for adidx = 1:2
            for fidx = 1:6
                
                nFolderIdx = (midx-1)*numel(asAngleDist)*numel(aAngleNum) + ...
                             anidx*numel(asAngleDist) + ... 
                             adidx;
                sFolderName = [asMode{midx}, '_Ntx', num2str(aAngleNum(anidx)), '_', asAngleDist{adidx},'\'];
                
                %%% Make folders
                if ~isdir(sFolderName)
                    mkdir(sFolderName);
                end
                
                %%% Gen Angles
                if strcmp(asAngleDist{adidx}, 'Equal_alpha')
                    stTxAngle.bEqualAlpha = 1;
                    stTxAngle.aAzi_deg = asind( linspace(sind(-51.9525), sind(51.9525), aAngleNum(anidx))' ); % equal alpha
                else 
                    stTxAngle.bEqualAlpha = 0;
                    stTxAngle.aAzi_deg = linspace(-51.9525, 51.9525, aAngleNum(anidx))'; % equal theta
                end
                
                save([sFolderName 'stTxAngle'],'stTxAngle');
                
            end
        end
    end
end

                       