% function [I_q,q,qerr_byq,qSort,I_phi,phivalues,pctExcitation,angleMap,phiMap] = radialAverage(experiment,runnum,images,L,D,x0,y0,dq,dphi,coors,A,goodPixels)
function [I_q,I_q_adu,q,qerr_byq,qSort,I_phi,phivalues,pctExcitation,angleMap,phiMap] = radialAverage(experiment,runnum,images,images_adu,L,D,x0,y0,dq,dphi,coors,A,goodPixels) %bms, changed input/output

%Constructs a weighted radial average of an image in q-space.
%           image :  array of photons or counts
%               L :  wavelength (in nm)
%               D :  distance from the Be window to the face of the detector
%           x0,y0 :  center of the detector relative to the map in 'coors'
%              dq :  width of q bands
%           coors :  array of pixel centers
%               A :  array of pixel areas
%         goodPixels :  an array of 1s & 0s indicating which pixels "count"
%         distMap :  for each pixel, the distance the scattering point is
%                    from the laser focus

% We don't use these forms b/c they don't take the physical locations of
% the pixels into account.

% qMap: array of pixel positions in q-space
% qMap = genQ(coors, L, D, x0, y0);

% oMap: array of subtended angles per pixel
% oMap = genAngle( coors, A, D, x0, y0 );

%% This section uses James' 'physicalQ' to find physically grounded values for path-length, lMap, position in q-space, q, and error in q, qerr.
MapRadii= @(x,y) sqrt(x.^2+y.^2);
radiiMap = MapRadii(coors.x-x0,coors.y-y0);
MapPhi = @(x,y) atand(y./x);
phiMap = MapPhi(coors.x-x0,coors.y-y0);phiMap(phiMap<0)=180+phiMap(phiMap<0);phiMap(coors.y-y0<0)=phiMap(coors.y-y0<0)+180;

%% PhysicalQ turns a list of physical coordinates into a list of q-based coordinates:
[lList,qList,qerr,diList,avg_angleList,distances] = cache_results(@physicalQ,{experiment,radiiMap(:),D,L});
% [lList,qList,qerr,diList,avg_angleList,distances] = physicalQ(experiment,radii(:),D,L);
lMap     =  reshape(lList,[388 185 32]);
qMap     =  reshape(qList,[388 185 32]);
diMap    =  reshape(diList,[388 185 32]);
angleMap =  reshape(avg_angleList,[388 185 32]);
focdistMap  =  reshape(distances.distfromfocus,[388 185 32]);
XeBeamLengthMap = reshape(distances.distintochamber,[388 185 32]);

%% Beryllium window and air attenuation corrections:
% For these, be sure to use T = I/I0 = e^-Epsilon*L(microns!)*C

if strcmp(experiment,'i0613') % Fix values before re-enabling.
%     Be_EpsilonC = 1.1650e-04; % Using the calculated value for 7100eV. http://henke.lbl.gov/optical_constants/filter2.html
%     Air_EpsilonC = 7.155e-06; % Using the calculated value for 7100eV. http://henke.lbl.gov/optical_constants/gastrn2.html
elseif strcmp(experiment,'b0114')
    Be_EpsilonC = 1.7083e-04;   % T=.91813 for 500microns. Using the calculated value for 8300eV. http://henke.lbl.gov/optical_constants/filter2.html
    Air_EpsilonC = 1.0121e-06;  % Using the calculated value for 8300eV. http://henke.lbl.gov/optical_constants/gastrn2.html
    Xe_EpsilonC = 2.1826e-06;   % 1cm of Xenon at 10.9 Torr transmits .97841 of 8300eV. This is the only gas that noticably attenuates the x-ray.
elseif strcmp(experiment,'56012') % Fix values before re-enabling.
%     Be_EpsilonC = 1.5351e-05; % Using the calculated value for 20100eV. http://henke.lbl.gov/optical_constants/filter2.html
%     Air_EpsilonC = 3.43e-07;  % Using the calculated value for 20100eV. http://henke.lbl.gov/optical_constants/gastrn2.html
end

BeLengthMap = 500/cosd(angleMap); % This is the length of Beryllium the photons pass through for each pixel.
attenuators.BeTransmissionMap = exp(-Be_EpsilonC*BeLengthMap);
AirLengthMap = D/cosd(angleMap); % This is the length of air the photons pass through for each pixel.
attenuators.AirTransmissionMap = exp(-Air_EpsilonC*AirLengthMap);

if strcmp(experiment,'b0114') && runnum==179
    XeScatterLengthMap = (diMap-D)./cosd(angleMap); % This is the length of Xenon the photons pass through AFTER scattering.
    XeScatterLengthMap(diList<(D+500))=0;
    attenuators.XeScatterTransmissionMap = exp(-Xe_EpsilonC*XeScatterLengthMap);
%     disp('Compensating for Xenon absorption after scattering.')
    
    % This is the length of Xenon the photons pass through BEFORE scattering.
    attenuators.XeBeamTransmissionMap = exp(-Xe_EpsilonC*XeBeamLengthMap);
%     disp('Compensating for Xenon absorption before scattering.')
else
    attenuators.XeScatterTransmissionMap = ones(size(AirLengthMap)); % Creates a dummy for non-xenon runs.
    attenuators.XeBeamTransmissionMap = ones(size(AirLengthMap)); % Creates a dummy for non-xenon runs.
    XeScatterLengthMap = ones(size(AirLengthMap));
end
% max_Air_transmission = max(attenuators.AirTransmissionMap(:))
% min_Air_transmission = min(attenuators.AirTransmissionMap(:))
% max_Be_transmission = max(attenuators.BeTransmissionMap(:))
% min_Be_transmission = min(attenuators.BeTransmissionMap(:))
% max_Xe_transmission = max(attenuators.XeTransmissionMap(:))
% min_Xe_transmission = min(attenuators.XeTransmissionMap(:))

%% Convert Area into effective Area, which accounts for the pixel not being perpendicular to the photon trajectory:
effectivelengths=abs(sind(90-angleMap)*110);
A = (A/110).*effectivelengths;

%% oMap: array of subtended angles per pixel
oMap = (A./(diMap.^2)).*(1-2*((L*qMap)/(4*pi)).^2).^3;
% oMap = (A/D^2).*(1-2*((L*qMap)/(4*pi)).^2).^3;

%% pixMap: cellarray of arrays of indicies of pixels
% q: array of centers of bands in q
% [qSort,q,phiSort,phivalues,polSorts] = binPix(qMap,dq,phiMap,dphi,goodPixels,angleMap); % base case
[qSort,q,phiSort,phivalues,polSorts] = cache_results(@binPix,{qMap,dq,phiMap,dphi,goodPixels,angleMap}); % caching version

%% Get qerr_bypixel and convert it to one value per q:
qerr_byq=zeros(size(q));
% for loop=5:length(q);
%     qerr_byq(loop)=max(qerr(pixMap{1,loop}));
% end
%
I_q=zeros([length(images(1,1,1,:)) length(q)]);
I_phi=zeros([length(images(1,1,1,:)) length(phivalues)]);

%% Adds all similar pixels to compute I:
for imagenumber=1:length(images(1,1,1,:));
%     [I_q(imagenumber,:),I_phi(imagenumber,:)] = constructIq(images(:,:,:,imagenumber),qMap,oMap,qSort,phiSort,polSorts,phiMap,lMap,diMap,angleMap,attenuators,A,radiiMap,q,phivalues); 
    [I_q(imagenumber,:),I_q_adu(imagenumber,:),I_phi(imagenumber,:)] = constructIq(images(:,:,:,imagenumber),images_adu(:,:,:,imagenumber),qMap,oMap,qSort,phiSort,polSorts,phiMap,lMap,diMap,angleMap,attenuators,A,radiiMap,q,phivalues); %bms, edited input/output 
end

%% Plot histograms for each q (temporary):
% figure(17);
% for thisq=1:length(q)
%     plot(pix_hist_x,pix_hists(thisq,:));legend(num2str(q(thisq)));
%     pause(1);
% end

%% For each q, find the average distance from the focus (for excitation percentage reasons.)
for i=1:length(q)
    distfromfocus(i)=mean(focdistMap(qSort{i}));
end

% figure(83);plot(q,distfromfocus);ylabel('Distance from focus, microns');xlabel('q');

%% Also find (for theoretical calculation) the distance to detector for each q
% for i=1:length(q)
%     angle_of_q(i)=mean(avg_angleList(qSort{i}));
%     temp_distintochamber(i)=mean(distances.distintochamber(qSort{i}));
%     temp_XenonBeamTransmission(i)=mean(attenuators.XeBeamTransmissionMap(qSort{i}));
%     temp_XeScatterpath(i)=mean(XeScatterLengthMap(qSort{i}));
%     temp_XeScatterTransmission(i)=mean(attenuators.XeScatterTransmissionMap(qSort{i}));
% end
% figure(107);plot(q(:),temp_distintochamber(:),q(:),temp_XeScatterpath(:));legend('Xe Path Before Scatter','Xe Path After Scatter');
% figure(108);plot(q(:),temp_XenonBeamTransmission(:),q(:),temp_XeScatterTransmission(:));legend('Xe Beam Transmission','Xe Scatter Transmission');
% [pctExcitation] = excitationprofile(distfromfocus);
[pctExcitation] = [];
% save('ExcitationInfo.mat','q','distfromfocus');

end

