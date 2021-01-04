% Sua Bae
%
% 2016-12-04
%
% plot plane wave angles in theta, sin(theta) domain
%


function Plot_Angles(bPlotTransmitAngles, stTxAngle)

    if(bPlotTransmitAngles)

        nTxAngleNum = numel(stTxAngle.aAzi_deg);
        aTxAngle    = stTxAngle.aAzi_deg;
        aTxAlpha    = sind(aTxAngle);    
        dTxAlpha    = aTxAlpha(2) - aTxAlpha(1);    

        %%%%%%%%%%%%%%%% plot
        figure('Position',[300, 100, 500, 300]);
        stem(aTxAngle, aTxAlpha, 'x-');
        for aidx = 1:nTxAngleNum
            hold on; line([0,aTxAngle(aidx)],[aTxAlpha(aidx),aTxAlpha(aidx)],'Color','k');
        end
        ylim([-0.8,0.8]);
        ylabel('alpha');
        xlabel('angle (deg)');
    end
    
end