% function [I_q,I_phi] = constructIq(image,qMap,oMap,qSort,phiSort,polSorts,phiMap,lMap,diMap,twothetaMap,attenuators,Area,radii,q,phivalues)
function [I_q,I_q_adu,I_phi] = constructIq(image,image_adu,qMap,oMap,qSort,phiSort,polSorts,phiMap,lMap,diMap,twothetaMap,attenuators,Area,radii,q,phivalues) %bms, changed input/output

%Converts pixel map to I(q), the intensity as a function of q.
%   image : array of photons or counts
%   oMap  : array of subtended angles per pixel
%   pixMap: cellarray of arrays of indicies of pixels
%   lMap  : array of exposed path lengths
%   Intensities : an array of intensities at each q value

% flatten image, angle map, and length map
iFlat = image(:);
iFlat_adu = image_adu(:); %bms, added line
oFlat = oMap(:);
phiFlat = phiMap(:);
lFlat = lMap(:);
R2Flat = diMap(:).^2 + radii(:).^2;
BeFlat = attenuators.BeTransmissionMap(:);
AirFlat = attenuators.AirTransmissionMap(:);
AreaFlat = Area(:);
twothetaFlat = twothetaMap(:);
XeFlat = attenuators.XeScatterTransmissionMap(:).*attenuators.XeBeamTransmissionMap(:);
polarizationFlat = (sind(phiFlat).^2) + (cosd(phiFlat).^2).*(cosd(twothetaFlat).^2);
polarizationFlat=ones(size(XeFlat));


%% anonymous function (thank goodness) that computes the intensity. Each pixel is normed by solid angle and path length exposed.
% radialSum = @(i) sum(iFlat(i).*oFlat(i).*lFlat(i))/sum(oFlat(i).*lFlat(i));
% radialSum = @(i) sum(iFlat(i).*oFlat(i))/sum(oFlat(i).*lFlat(i)); %current best
% radialSum = @(i) sum(iFlat(i).*oFlat(i).*(1./lFlat(i)))/sum(oFlat(i).*(1./lFlat(i)));
% radialSum = @(i) sum(iFlat(i))/sum(oFlat(i).*lFlat(i)); %Vale's way

% radialSum = @(i) mean(iFlat(i)./(oFlat(i).*lFlat(i).*BeFlat(i).*AirFlat(i))); % Jim's way
pixelAverage = @(i) mean((iFlat(i).*R2Flat(i))./(AreaFlat(i).*lFlat(i).*BeFlat(i).*AirFlat(i).*XeFlat(i).*polarizationFlat(i))); % Jim's way (Xenon is ones for nonxenon runs).
pixelSum = @(i) sum(iFlat_adu(i)); %bms, added line. sum of all ADU in specified q range
% radialSum = @(i) mean(iFlat(i)./(Area(i).*BeFlat(i).*AirFlat(i).*XeFlat(i))); % What does it look like without path length correction?
% radialSum = @(i) mean(iFlat(i)./(lFlat(i).*BeFlat(i).*AirFlat(i).*XeFlat(i))); % How about without dividing by area?
% radialSum = @(i) mean(iFlat(i)./(oFlat(i).*BeFlat(i).*AirFlat(i))); % Alternate testing method, which ignores reaction length
% pixelSum  = @(i) mean(iFlat(i));
% radialSum = @(i) sum(iFlat(i)./lFlat(i))/sum(oFlat(i));

I_q=cellfun(pixelAverage,qSort);
I_q_adu=cellfun(pixelSum,qSort); %bms, added line
I_phi=cellfun(pixelAverage,phiSort); % Does the same process, just grouped by phi now instead of q.

% I_0deg=cellfun(pixelAverage,polSorts.inplane);
% I_45deg=cellfun(pixelAverage,polSorts.unpolarized);
% I_90deg=cellfun(pixelAverage,polSorts.perpendicular);
% figure(90);
% plot(q,I_0deg,q,I_45deg,q,I_90deg);
% legend('In-Plane with Polarization','45 degrees from polarization','Perpendicular to polarization');xlim([1 4.5]);

I_1A=cellfun(pixelAverage,polSorts.oneA);
I_2A=cellfun(pixelAverage,polSorts.twoA);
I_3A=cellfun(pixelAverage,polSorts.threeA);
I_4A=cellfun(pixelAverage,polSorts.fourA);
figure(91);
subplot(2,2,1);imagesc(CsPadRearrangeXPP(reshape(polarizationFlat,[388,185,32])));title('Polarization Factor');axis square;colorbar;
subplot(2,2,3);imagesc(CsPadRearrangeXPP(reshape(phiFlat,[388,185,32])));title('Phi');axis square;colorbar;
subplot(1,2,2);plot(phivalues,I_2A,phivalues,I_3A,phivalues,I_4A);
legend('2 Angstrom','3 Angstrom','4 Angstrom');xlim([0 360]);

for j=1:length(qSort)
    try
        mean_i(j)=mean(iFlat(qSort{j}));
        mean_di(j)=mean(diMap(qSort{j}));
        mean_radius(j)=mean(radii(qSort{j}));
        mean_R2(j)=mean(R2Flat(qSort{j}));
        mean_area(j)=mean(Area(qSort{j}));
        mean_l(j)=mean(lFlat(qSort{j}));
        mean_be(j)=mean(BeFlat(qSort{j}));
        mean_air(j)=mean(AirFlat(qSort{j}));
        mean_xe(j)=mean(XeFlat(qSort{j}));
        mean_polarization(j)=mean(polarizationFlat(qSort{j}));
        mean_angleFlat(j)=mean(twothetaFlat(qSort{j}));
    catch
        mean_i(j)=0;
        mean_di(j)=0;
        mean_radius(j)=0;
        mean_R2(j)=0;
        mean_area(j)=0;
        mean_l(j)=0;
        mean_be(j)=0;
        mean_air(j)=0;
        mean_xe(j)=0;
        mean_polarization(j)=0;
        mean_angleFlat(j)=0;
    end
end

mean_i=mean_i;
mean_di=mean_di/max(mean_di);
mean_radius=mean_radius/max(mean_radius);
mean_R2=mean_R2/max(mean_R2);
mean_area=mean_area/max(mean_area);
mean_l_normalized=mean_l/max(mean_l);
mean_be=mean_be/max(mean_be);
mean_air=mean_air/max(mean_air);
mean_xe=mean_xe/max(mean_xe);
mean_polarization=mean_polarization/max(mean_polarization);
mean_angleFlat=mean_angleFlat/max(mean_angleFlat);

figure(852);
plot(q,mean_R2,q,mean_area,q,mean_l_normalized,q,mean_be,q,mean_air,q,mean_xe,q,mean_polarization,q,mean_angleFlat,q,mean_di,q,mean_radius);
legend('R^2','Effective pixel area','Path length','Beryllium transmission','Air Transmission','Xenon Transmission','Polarization factor','Angle','di','radius');
save('allCorrectionFactors.mat','mean_i','mean_R2','mean_area','mean_l','mean_be','mean_air','mean_xe','mean_polarization','mean_angleFlat','mean_di','q','mean_radius');

save('pixMap2_temp.mat','qMap','qSort','q','phiMap','phiSort','phivalues','polSorts','phiFlat','twothetaMap','oFlat','lFlat','R2Flat','BeFlat','AirFlat','AreaFlat','XeFlat');

end

