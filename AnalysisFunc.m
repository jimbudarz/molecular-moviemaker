function [binned_img,binnedN,binnedIPM,photon_img,sum_img,N_shots_out,in_struc]=AnalysisFunc(in_struc)
%Applies gain and common mode corrections and filters, and then creates a
%photon map.

%% Rename appropriate variables
currentshots=in_struc.currentshots;
currentstep=in_struc.step;
isXRayon=in_struc.isXRayon;
isphotoncounting=in_struc.isphotoncounting; % 0 is not photon counting, 1 is photon counting, .5 is the new hybrid method.
if isXRayon==1
    fina=in_struc.fina;
elseif isXRayon==0
    fina=in_struc.darkfina;
end

%% Read in CSPAD data
% [df,CSPADtsMS]=rdCSPADdataXPP(fina,currentstep,currentshots);
% CSPADtsMS = CSPADtsMS(currentshots);
% % Note that the CSPAD timestamps are ALL pulled at first (only reading
% % a subset causes a mismatch because it subtracts the minimum one. 

%% Read in CSPAD data in parallel!
% First divide shots by into subsets:
if matlabpool('size')==0
    [df,CSPADtsMS]=rdCSPADdataXPP(fina,currentstep,currentshots);
else
    [df,CSPADtsMS]=rdCSPADdataParallel(fina,currentstep,currentshots);
end
CSPADtsMS = CSPADtsMS(currentshots);
% Note that the CSPAD timestamps are ALL pulled at first (only reading
% a subset causes a mismatch because it subtracts the minimum one. 

%% Use the dark shots to make a pedestal and a pedestal mask

if isXRayon==1 % Subtract pedestal
    df=df-repmat(in_struc.pedestal,[1,1,1,size(df,4)]);
end

%% Apply gain correction:
if strcmp(in_struc.experiment,'56012')
    df=df.*in_struc.gain_corr(:,:,:,1:size(df,4));
end

%% Common Mode Correction, per-ASIC basis (only works if most pixels are not lit up):

if strcmp(in_struc.experiment,'56012') && isXRayon==1
    cModeMask=in_struc.pedestalmask;
    for shot = 1:numel(df(1,1,1,:))
        [df(:,:,:,shot),~]=corrCNoise(df(:,:,:,shot),cModeMask,'budarz');
    end
elseif (strcmp(in_struc.experiment,'b0114') || strcmp(in_struc.experiment,'i0613')) && isXRayon==1
    cModeMask=in_struc.pedestalmask.*in_struc.noisemask;
    for shot = 1:size(df,4)
        [df(:,:,:,shot),~]=corrCNoise(df(:,:,:,shot),cModeMask,in_struc.sigma,'sigma');
%         [df(:,:,:,shot),~]=corrCNoise(df(:,:,:,shot),cModeMask,in_struc.sigma,'budarz');
    end
end

%% Now that common mode is done, apply all filters:
if isXRayon==1
    df = df.*in_struc.repGoodPixels(:,:,:,1:size(df,4));
end

%% Take straight sum of images without photon counting:
sum_img=sum(df,4);

%% Set up time-tool data for the current set of shots:

if isXRayon==1
    if in_struc.isUVon == 1;
        currentTTdata_sec=in_struc.IPMandTTdata.goodTTdata_sec(in_struc.minshot:in_struc.maxshot);
        currentTTts=in_struc.IPMandTTdata.goodTTtsMS(in_struc.minshot:in_struc.maxshot);
        in_struc.currentRealTime_sec=in_struc.current_timestep_sec+currentTTdata_sec;
%         figure(253);subplot(1,2,1);hist(in_struc.currentRealTime_sec,50);subplot(1,2,2);plot(in_struc.currentRealTime_sec);pause(1);
    elseif in_struc.isUVon == 0;
        currentTTdata_sec=in_struc.IPMandTTdata.noUVTTdata_sec(in_struc.minshot:in_struc.maxshot);
        currentTTts=in_struc.IPMandTTdata.noUVTTtsMS(in_struc.minshot:in_struc.maxshot);
        in_struc.currentRealTime_sec=repmat(in_struc.current_timestep_sec,size(currentTTdata_sec)); % Sets time-tool offset for UV-off frames to 0.
    end
    % Verify that the right time-tool data is applied to the right frame:
    try
        if currentTTts ~= CSPADtsMS
            disp('Time-tool timestamps do not match CSPAD timestamps.');
        end
    catch
        disp('Time-tool data does not match shots! Check TTdataERROR.mat.');
        TTtsMS=in_struc.IPMandTTdata.TTtsMS;
        save('TTdataERROR.mat','currentshots','CSPADtsMS','currentTTts','TTtsMS');
    end
end

%% Make a map of photons

if isphotoncounting==1 && not(isempty(df))
    [counts_img,photon_img,N_shots_out]=photon_map(df,in_struc);
    binned_img=0;binnedN=0;in_struc.binnedcountavg=0;in_struc.realTimes=0;
elseif isphotoncounting==0
    counts_img=sum_img;photon_img=sum_img;N_shots_out=size(df,4);
    binned_img=0;binnedN=0;in_struc.binnedcountavg=0;in_struc.realTimes=0;
elseif isphotoncounting==.5
    if isXRayon==1
        [binned_img,binnedN,binnedIPM,photon_img,N_shots_out,in_struc.count_sums,in_struc.binnedcounttot,in_struc.realTimes]=photon_map_hybrid(df,in_struc);
    elseif isXRayon==0
        photon_img=sum_img;N_shots_out=size(df,4); % The lower threshold of hybrid photon counting here would not make a good new pedestal.
        binned_img=0;binnedN=0;binnedIPM=0;
    end
end


end