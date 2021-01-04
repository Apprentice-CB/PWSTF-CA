%
% 2015-09-10
% Sua Bae
%
% Generate TX or RX source positions
%
% mSrcPos           : source poisition and element index, dim = (xpos, ypos, zpos, eleidx)x(nRxAptrEleNum*nSrcNumPerTransEle)
%
% sTransType        : 'Linear', 'Convex', 'Matrix'
% nTransElePitch    : [meter] element pitch,  [pitch_x, pitch_y] when 'Matrix' 
% nTransEleWidth    : [meter] element width,  [width_x, width_y] when 'Matrix' 
% nTransRadius      : [meter] radius of Xdcr when 'Convex'
% nAptrEleNum       : number of elements, [num_x, num_y] when 'Matrix'
% nSrcNumPerEle     : number of source per element, [num_x, num_y] when 'Matrix'
% nAptrCenterPos    : Position of center of source array
% nAptrCenterAngle  : Angle of center of source array
function mSrcPos = GenerateSource(sTransType, nTransElePitch, nTransEleWidth, nTransRadius, ...
                                    nTransEleNum, nSrcNumPerEle, nAptrCenterPos, nAptrCenterAngle)

    switch sTransType
        case 'Linear'
            aElePos    = -nTransElePitch*(nTransEleNum-1)/2:nTransElePitch:nTransElePitch*(nTransEleNum-1)/2; % x positions of elements
            aElePos    = aElePos + nAptrCenterPos;

            nSrcInterval = nTransEleWidth/nSrcNumPerEle;
%             nSrcInterval = nTransEleWidth/(nSrcNumPerEle-1); % -- wrong

            aSrcRelativePos = ( (-nSrcInterval*(nSrcNumPerEle-1)/2) : nSrcInterval : (nSrcInterval*(nSrcNumPerEle-1)/2) )'; % distance of source points from element position

            mSrcPos = zeros(4,nTransEleNum*nSrcNumPerEle);    
            for eidx = 1:nTransEleNum
                for sidx = 1:nSrcNumPerEle
                    mSrcPos(1,nSrcNumPerEle*(eidx-1)+sidx) = aElePos(eidx) + aSrcRelativePos(sidx); % x position
                    mSrcPos(2,nSrcNumPerEle*(eidx-1)+sidx) = 0;  % y position
                    mSrcPos(3,nSrcNumPerEle*(eidx-1)+sidx) = 0;  % z position
                    mSrcPos(4,nSrcNumPerEle*(eidx-1)+sidx) = eidx;  % index of elements which the source belongs to
                end
            end
            
        case 'Convex'
            aElePos_theta  = ( -nTransElePitch*(nTransEleNum-1)/2:nTransElePitch:nTransElePitch*(nTransEleNum-1)/2 )' / nTransRadius; % [rad] theta positions of elements            
%             figure;compass(exp(1i*aElePos_theta));
%             axis1 = gca;
%             axis1.CameraUpVector = [-1 0 0];
            
            nSrcInterval = nTransEleWidth/nSrcNumPerEle; % [meter]
            aSrcRelativeTheta = ( (-nSrcInterval*(nSrcNumPerEle-1)/2) : nSrcInterval : (nSrcInterval*(nSrcNumPerEle-1)/2) )' / nTransRadius; % [rad]  distance of source points from element position
           
            mSrcPos = zeros(4,nTransEleNum*nSrcNumPerEle);           
            for eidx = 1:nTransEleNum
                for sidx = 1:nSrcNumPerEle
                    nTheta = aElePos_theta(eidx) + aSrcRelativeTheta(sidx); % [rad]
                    mSrcPos(1,nSrcNumPerEle*(eidx-1)+sidx) = nTransRadius*sin(nTheta); % x position
                    mSrcPos(2,nSrcNumPerEle*(eidx-1)+sidx) = 0;  % y position
                    mSrcPos(3,nSrcNumPerEle*(eidx-1)+sidx) = nTransRadius*cos(nTheta);  % z position
                    mSrcPos(4,nSrcNumPerEle*(eidx-1)+sidx) = eidx;  % index of elements which the source belongs to
                end
            end
            
            figure;
            scatter(mSrcPos(1,:),mSrcPos(3,:));
            axis equal; axis tight; set(gca,'Ydir','reverse');ylim([0, nTransRadius + 120e-3]);xlim([-100e-3, 100e-3]); grid on;
            
        case 'Matrix'
            
            aElePos_x    = -nTransElePitch(1)*(nTransEleNum(1)-1)/2:nTransElePitch(1):nTransElePitch(1)*(nTransEleNum(1)-1)/2; % x positions of elements
            aElePos_y    = -nTransElePitch(2)*(nTransEleNum(2)-1)/2:nTransElePitch(2):nTransElePitch(2)*(nTransEleNum(2)-1)/2; % y positions of elements

            nSrcInterval_x = nTransEleWidth(1)/nSrcNumPerEle(1);
            nSrcInterval_y = nTransEleWidth(2)/nSrcNumPerEle(2);
%             nSrcInterval = nTransEleWidth/(nSrcNumPerEle-1); % -- wrong

            aSrcRelativePos_x = ( (-nSrcInterval_x*(nSrcNumPerEle(1)-1)/2) : nSrcInterval_x : (nSrcInterval_x*(nSrcNumPerEle(1)-1)/2) )'; % x distance of source points from element position
            aSrcRelativePos_y = ( (-nSrcInterval_y*(nSrcNumPerEle(2)-1)/2) : nSrcInterval_y : (nSrcInterval_y*(nSrcNumPerEle(2)-1)/2) )'; % y distance of source points from element position

            mSrcPos = zeros(5,nTransEleNum(1)*nTransEleNum(2)*nSrcNumPerEle(1)*nSrcNumPerEle(2));    

            for eidx_y = 1:nTransEleNum(2)            
                for eidx_x = 1:nTransEleNum(1)
                    for sidx_y = 1:nSrcNumPerEle(2)
                        for sidx_x = 1:nSrcNumPerEle(1)                        
                            idx =   (eidx_y-1)*nSrcNumPerEle(1)*nSrcNumPerEle(2)*nTransEleNum(1) ...
                                  + (eidx_x-1)*nSrcNumPerEle(1)*nSrcNumPerEle(2) ...
                                  + (sidx_y-1)*nSrcNumPerEle(1) ...
                                  + sidx_x;
                            mSrcPos(1,idx) = aElePos_x(eidx_x) + aSrcRelativePos_x(sidx_x); % x position
                            mSrcPos(2,idx) = aElePos_y(eidx_y) + aSrcRelativePos_y(sidx_y); % y position
                            mSrcPos(3,idx) = 0;  % z position
                            mSrcPos(4,idx) = (eidx_y-1)*nTransEleNum(1) + eidx_x; % index of elements which the source belongs to
                            mSrcPos(5,idx) = eidx_x;  % x index of elements which the source belongs to
                            mSrcPos(6,idx) = eidx_y;  % y index of elements which the source belongs to
                            mSrcPos(7,idx) = sidx_x;  % sub x index within a element
                            mSrcPos(8,idx) = sidx_y;  % sub y index within a element
                        end
                    end
                end
            end
    end
    
%     
%     mElePos = repmat(aElePos,nSrcNumPerEle,,
%     
%     element_center_point = - nAptrSize/2 +zeta : nTransElePitch : nAptrSize/2 + zeta;
%     nL_center_point = length(element_center_point);
%     
%     element_width_point = linspace(-nTransEleWidth/2, nTransEleWidth/2, width_sample_num);
%     
%     rep_element_point = repmat(element_center_point,[width_sample_num 1]);
%     rep_element_width_calc = repmat(element_width_point',[1 nL_center_point]);
%     
%     calc_element_point = rep_element_point + rep_element_width_calc;
%     aperture = reshape(calc_element_point, 1, nL_center_point*width_sample_num);
end