function [binned_img,binnedN,binnedIPM,photon_img,N_shots_out,count_sums,binnedcounttot,currentRealTime_sec]=photon_map_hybrid(df,in_struc)
%Uses a hybrid photon counting method (minimum ADU threshold) and then bins
%data according to time tool.

%% Initial Parameters
currentRealTime_sec=in_struc.currentRealTime_sec;
timebins_sec=in_struc.timebins_sec;
min_count=in_struc.min_count(:,:,:,1:size(df,4)); % Sets the lower boundary on the detector (in ADU) to define a photon event.
% gain_corr=repmat(in_struc.gain_corr,[1 1 1 shots_this_loop]);
% adu2photons = 24; %the average number of ADU that signal one photon

%% Apply gain correction factors (Currently disabled because good gain correction for i0613 hasn't been made.)
% df=df.*gain_corr.*isfinite(gain_corr); % df is a matrix of the detector images that have been read in.  The first 3 dimensions define the pixel, the 4th defines the frame number.

%% Make a plot of photons
% myphotons=df.*repmat(in_struc.goodMap,[1,1,1,size(df,4)]); % Because of count_sums, only useful pixels should be accounted for after this point.
myphotons=df; % The above line makes no difference. Maybe AnalysisFunc takes care of it.

% min_count=-100.*ones(size(min_count));
myphotons((df>-min_count & df<min_count) | not(isfinite(myphotons)))=0; 
% myphotons=myphotons.*(df>min_count & isfinite(myphotons)); % Uses only a single lower limit to define a photon or partial photon (above the 0-photon peak)
% myphotons(df<(min_count) & df>(-1*min_count))=0; % This method would allow the counting of very low values from random fluctuation:

%% For each frame, sum over the whole detector to estimate X-Ray flux:
count_sums = sum(sum(sum(myphotons,3),2),1);
count_sums = count_sums(:);

%% Plot count histogram
% if currentRealTime_sec==zeros(size(currentRealTime_sec))
%     figure(201);subplot(1,1,1);hist(count_sums,50);title('Count Sums');
% else
%     figure(201);subplot(1,2,1);hist(count_sums,50);title('Count Sums');subplot(1,2,2);hist(currentRealTime_sec,50);title('TT Jitter');pause(1);
% end

%% (Optional) Discard frames with abnormally high or low total intensities:
% tossableshots=count_sums<(mean(count_sums)*.8) | count_sums>(mean(count_sums)*.9);
% % tossableshots=(count_sums<in_struc.count_sum_low_limit | count_sums>in_struc.count_sum_high_limit);
% disp([num2str(100*sum(tossableshots)/size(myphotons,4)) '% of shots thrown away for high or low intensity.']);
% 
% count_sums(tossableshots) = [];
% myphotons(:,:,:,tossableshots) = [];
% currentRealTime_sec(tossableshots) = [];
% disp(['      Average of ' num2str(mean(count_sums),'%.2e') ' total ADU per frame.']);

%% (Optional) Scale the frames by their count totals: (MOVED TO LARGERUN AFTER BINNING DUE TO NOISE/STATISTICS CONCERNS)
% myphotons = myphotons./permute(repmat(count_sums,[1 388 185 32]),[2 3 4 1]);

%% Sum all images (unsorted by time-tool)
photon_img = sum(myphotons,4);
% pixel_img_photons = photon_img./adu2photons; %sum of all images (unbinned) in photons

%% Sort images according to time tool data:

% Preallocate
binned_img=zeros([388 185 32 length(timebins_sec)]);   % Binned data will now have the 4th dimension as the time-bin.
% binned_pixel_img_photons=zeros([388 185 32 length(timebins_sec)]);   % Binned data will now have the 4th dimension as the time-bin.
binnedN=zeros([length(timebins_sec) 1]);               % Binned data, as a sum, requires a corresponding number of images so an average can be taken later.
binnedcounttot=zeros([length(timebins_sec) 1]);      % Binned detector sums for troubleshooting.
binnedIPM=zeros([length(timebins_sec) 1]);

% Bins are defined as extending halfway to the next nominal bin position:    
for binnum=1
    binmin = timebins_sec(binnum);
    binmax = (timebins_sec(binnum+1)+timebins_sec(binnum))/2;
    matchingshots = find(currentRealTime_sec<binmax & currentRealTime_sec>binmin);
    binned_img(:,:,:,binnum) = sum(myphotons(:,:,:,matchingshots),4);
%     binned_pixel_img_photons(:,:,:,binnum) = sum(myphotons(:,:,:,matchingshots),4)./adu2photons;
    binnedN(binnum) = length(matchingshots); %number of matching shots per bin
    binnedcounttot(binnum)=sum(count_sums(matchingshots));
    binnedIPM(binnum) = sum(in_struc.IPMandTTdata.IPMdata(matchingshots));
    clear binmin binmax matchingshots;
end

for binnum=length(timebins_sec)
    binmin = (timebins_sec(binnum)+timebins_sec(binnum-1))/2;
    binmax = timebins_sec(binnum);
    matchingshots = find(currentRealTime_sec<binmax & currentRealTime_sec>binmin);
    binned_img(:,:,:,binnum) = sum(myphotons(:,:,:,matchingshots),4);
%     binned_pixel_img_photons(:,:,:,binnum) = sum(myphotons(:,:,:,matchingshots),4)./adu2photons;
    binnedN(binnum) = length(matchingshots);
    binnedcounttot(binnum)=sum(count_sums(matchingshots));
    binnedIPM(binnum) = sum(in_struc.IPMandTTdata.IPMdata(matchingshots));
    clear binmin binmax matchingshots;
end

for binnum=2:length(timebins_sec)-1
    binmin(binnum) = (timebins_sec(binnum)+timebins_sec(binnum-1))/2;
    binmax(binnum) = (timebins_sec(binnum+1)+timebins_sec(binnum))/2;
    matchingshots{binnum} = find(currentRealTime_sec<binmax(binnum) & currentRealTime_sec>binmin(binnum));
    binned_img(:,:,:,binnum) = sum(myphotons(:,:,:,matchingshots{binnum}),4);
%     binned_pixel_img_photons(:,:,:,binnum) = sum(myphotons(:,:,:,matchingshots),4)./adu2photons;
    binnedN(binnum) = length(matchingshots{binnum});
    binnedcounttot(binnum)=sum(count_sums(matchingshots{binnum}));
    binnedIPM(binnum) = sum(in_struc.IPMandTTdata.IPMdata(matchingshots{binnum}));
end

%% Keep track of the total number of images that have been summed (for unsorted data)
N_shots_out=length(count_sums);

end
