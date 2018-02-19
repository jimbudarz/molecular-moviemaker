function [I_allsteps_good_unbinned,q,photonimg_good_unbinned_avg,photonimg_noUV_unbinned_avg,IPMandTTdata]=LargeRun(experiment,runnum,timebins_fs,N_shots,steps)
% This function processes the HDF5 files from pump-probe LCLS experiments which employ the CS-PAD imaging detector.
% It is particularly designed to produce radial averages of 2D images as a
% function of time with all necessary calibrations and adjustments: I(q,t) 
%
%    experiment : The proposal number for the experiment (e.g. 'L560','560','l560','56012')
%        runnum : The run (a.k.a. scan) number to be loaded and processed 
%   timebins_fs : An array of the desired timebins to produce from the timepoints in the scan (uses time-tool to correct UV/Xray jitter) 
%       N_shots : The number of shots to load. Defaults to all available.
%         steps : The number of timepoints to load. Defaults to all available.
%
%      I_allsteps_good_unbinned : I(q,t), Intensity for laser-on shots (fordiagnostics)
%                             q : Array of momentum-transfer vector values corresponding to I(q,t)
%   photonimg_good_unbinned_avg : Average image for all laser-on frames
%   photonimg_noUV_unbinned_avg : Average image for all laser-off frames
%                  IPMandTTdata : Intensity+Position monitor data and time-tool (jitter correction)
%
% Note: the most important results generated here are saved in .mat files
% for further processing and are not necessarily the output of this
% function.

%% Preliminary setup of variables
try
    matlabpool open
catch
    disp([num2str(matlabpool('size')) ' workers already open.']);
end
runstarttime=tic; my_time_stamp=datestr(now,'yyyy-mmm-dd-HHMMSS'); disp(['Starting Run ',num2str(runnum)]);
shots_per_loop=50;  % To prevent overloading memory, frames are loaded in batches.
in_struc.count_sums_all=[];in_struc.ipm_sums_all=[];countsums_good_unbinned=[];countsums_noUV_unbinned=[];countsums_binned_tot=0;realTimes=[];

%% Set any experiment-specific variables:
switch experiment
    case {'L560','560','l560','56012'}
        experiment='56012'; hutch='xpp';
        in_struc.min_count = 110*ones([388 185 32 shots_per_loop]);
        STDlimit=3;
        fina=GetExpFina(hutch,experiment,runnum);
        D=40123; x0 = 94370; y0 = 92520;% Distance between detector and beryllium window holder.
        in_struc.isphotoncounting=.5;
        plottingqmin=2; plottingqmax=10;
        load('L560_masks.mat');%in_struc.pedestal=repmat(pedestal,[1,1,1,shots_per_loop]);
        load('gain_corr_L560.mat');in_struc.gain_corr = repmat(gain_corr,[1,1,1,shots_per_loop]);clear gain_corr;
        goodPixels=pedestalmask.*hotmask.*deadmask; %goodPixels=ones([388 185 32]);
        in_struc.goodPixels=goodPixels; in_struc.repGoodPixels = repmat(in_struc.goodPixels,[1,1,1,shots_per_loop]);in_struc.hotmask=hotmask;
        in_struc.pedestalmask=pedestalmask;
        timebins_sec = (1e-15).*[-20000:1000:5000];
        
    case {'i0613','xppi0613'}
        experiment='i0613'; hutch='xpp';
        fina=GetExpFina(hutch,experiment,runnum);
        plottingqmin=.75; plottingqmax=3.65; % For plotting, set the boundaries we expect to be accurate.
        load('i0613_darkstats');
        STDlimit=3;
        in_struc.min_count = 2*allsigmas; in_struc.min_count=repmat(in_struc.min_count,[1,1,1,shots_per_loop]);
        D=35972; x0 = 94370; y0 = 92520; % Measured distance between detector and beryllium window holder.
        load('i0613_masks.mat'); in_struc.pedestal=repmat(pedestal,[1,1,1,shots_per_loop]); % Loads masks of bad pixels and pedestal.
        in_struc.isphotoncounting=.5; % .5 refers to the hybrid photon counting method that only relies on a lower limit for "hits".
        in_struc.hotmask=hotmask;in_struc.noisemask=noisemask_strict;
        goodPixels=pedestalmask.*noisemask_lax.*hotmask.*deadmask; in_struc.goodPixels=goodPixels; in_struc.repGoodPixels = repmat(in_struc.goodPixels,[1,1,1,shots_per_loop]);
        if ismember(runnum,258:266)
            in_struc.count_sum_low_limit=3e7;
            in_struc.count_sum_high_limit=5e7;
            timebins_sec = (1e-15).*[-10200;-10100;-10000;-9900;-9800;-5200;-5100;-5000;-4900;-4800;-2200;-2100;-2000;-1900;-1800;-500;-400;-300;-200;-100;0;100;200;300;400;500;1800;1900;2000;2100;2200;4800;4900;5000;5100;5200;9800;9900;10000;10100;10200];
        elseif ismember(runnum,[225:257,267:325])
            in_struc.count_sum_low_limit=1.5e7;
            in_struc.count_sum_high_limit=2.5e7;
            timebins_sec = (1e-15).*[-2000:100:2000];
        elseif ismember(runnum,210)
            in_struc.count_sum_low_limit=1.3e7;
            in_struc.count_sum_high_limit=1.5e7;
            timebins_sec = (1e-15).*[-10000;-5000;-2500;-1000;-800;-600;-500;-400;-300;-200;-100;0;100;200;300;400;500;600;700;800;900;1000;1500;2000;3000;5000;10000];
        end
        
    case {'LB01','lb01','b01','b0114'}
        experiment='b0114';hutch='xpp';
        fina=GetExpFina(hutch,experiment,runnum);
        if runnum<258
            D=38770;
        elseif runnum>257
            D=38200;
        end        
        x0 = 95020; y0 = 92670; % Distance between detector and beryllium window holder and calculated image centers (not cspad center)
        plottingqmin=1.0; plottingqmax=4.2; % For plotting, set the boundaries we expect to be accurate.
        STDlimit=3;
        load('b0114darkmasks.mat'); in_struc.min_count = repmat(2*allsigmas,[1 1 1 shots_per_loop]); in_struc.sigma=allsigmas;
        in_struc.isphotoncounting=.5;
        goodPixels=noisemask.*pedestalmask;in_struc.goodPixels=goodPixels;in_struc.repGoodPixels=repmat(in_struc.goodPixels,[1,1,1,shots_per_loop]);
        pedestalmask=goodPixels;
        in_struc.noisemask=noisemask;in_struc.hotmask=ones([388,185,32]);
        in_struc.count_sum_low_limit=0;
        in_struc.count_sum_high_limit=Inf;
        timebins_sec = (1e-15).*timebins_fs;
        if ismember(runnum,[40,177,258,260,262])
            isdarkrun=1;
        else
            isdarkrun=0;
        end
        
    otherwise
        disp('Failed to recognize experiment');
end

%% Fetch run-specific parameters
[timesteps_sec]=fetch_timepoints(fina);     % Loads all preset timepoints for current run.
% timesteps_sec = -1*timesteps_sec;          % Fixes XPP's convention that negative time were UV-early (opposite of our convention).
[ordered_timepoints_sec,timeorder_steps]=sort(timesteps_sec,'ascend'); % Sorts random timepoints into their proper order.
totalnumsteps=length(timesteps_sec);        % Sets the total number of steps found.
[~,cspadtsinms]=rdCSPADTS(fina,1); totalframesperstep=length(cspadtsinms); % Sets the total number of shots per step based on the number in the first step.
allpressures=zeros([totalnumsteps 1]);      % Allocates a vector to later contain the average pressure during each step.

% If the desired step is not defined as an input, all steps are used.
if nargin<5
    steps=1:totalnumsteps;
end

% If the desired number of frames is not defined as an input, all good frames are used.
if nargin<4
    N_shots=totalframesperstep;
end

%% Rename some values to be passed into AnalysisFunc
in_struc.fina=fina;in_struc.darkfina=fina;
in_struc.shots_per_loop=shots_per_loop;
in_struc.experiment=experiment;
in_struc.timebins_sec=timebins_sec;
in_struc.D=D;in_struc.x0=x0;in_struc.y0=y0;
in_struc.goodPixels=goodPixels;
darksearch=0;

%% Parse shots into lists of only the desirable ones:
stepstoignore=zeros([length(totalnumsteps) 1]);
parfor steploop=1:totalnumsteps
    [goodshots{steploop},darkshots{steploop},noUVshots{steploop},IPMandTTdata{steploop}] = goodshotsonly(hutch,experiment,runnum,steploop,N_shots); % Parses the list of all shots into different types of shots.
end

% If there's no dark shots in this run, look for some in other runs.
while isempty(darkshots{1})
    darksearch=darksearch+1;disp(['Checking run ',num2str(runnum+darksearch),' for darks.'])
    in_struc.darkfina=GetExpFina(hutch,experiment,runnum+darksearch);
    parfor steploop=1:totalnumsteps
        [~,darkshots{steploop},~,~] = goodshotsonly(experiment,runnum+darksearch,steploop,N_shots);
    end
end;
disp(['Using darks from run ',num2str(runnum+darksearch)])

if isdarkrun==1
    darkshots{1}=1:N_shots;
else
    steps(stepstoignore==1)=[];
end

%% Fetch X-Ray wavelength for precise calculation of q:
Wavelength=zeros([1 length(steps)]);
for steploop=steps
    try
        [Wavelength(steploop),~] = getXRayWavelength(fina,steploop); % (Assigns lambda in nm.)
        disp(['X-Ray recorded as ' num2str(1.23984187/Wavelength(steploop)) ' keV.']);
    catch
        if strcmp(experiment,'56012')
            Wavelength(steploop)= 1.23984187/(20.1);
            disp(['Error: X-ray wavelength not found. Default of ',num2str(1.23984187/Wavelength(steploop)),' keV will be used.']);
        else
            disp('Error: X-ray wavelength not found.  Possible data damage.');
        end
    end
end
in_struc.Wavelength=Wavelength;

%% Make Pedestal from ALL steps of this run:
photonimg_dark_tot=zeros([388,185,32,length(totalnumsteps)]);N_darkshots_counter=zeros([length(totalnumsteps) 1]);in_struc.isXRayon=0;
% darkstarttime=tic;

for s=1:totalnumsteps
    if length(darkshots{s})>N_shots
        darkshots{s}=darkshots{s}(1:N_shots);
    end
    [photonimg_dark_tot(:,:,:,s),N_darkshots_counter(s)] = generatePedestal(experiment,darkshots{s},shots_per_loop,in_struc,s);
    % Too much input data to cache results:
    % [photonimg_dark_tot(:,:,:,s),N_darkshots_counter(s)] = cache_results(@generatePedestal,{experiment,currentdarkshots,shots_per_loop,in_struc,s});
end

% Create a pedestal from the dark image in ADU/shot.
photonimg_dark_avg=sum(photonimg_dark_tot,4)/sum(N_darkshots_counter);
 
if not(strcmp(experiment,'56012'))% || not(strcmp(runnum,'238'))
    pedestal=photonimg_dark_avg;
end

if N_darkshots_counter==0 % Find dark data when it's missing from within run.
    disp('No dark shots found, using most recent dark run instead!');
    if runnum>260
        pedestal=load(['~/cache/tempPedestals/',experiment,'/Run',num2str(260)],'pedestal');    
    else
        pedestal=load(['~/cache/tempPedestals/',experiment,'/Run',num2str()],'pedestal');
    end
else
    tempPedestalName=['~/cache/tempPedestals/',experiment,'/Run',num2str(runnum)];
    save(tempPedestalName,'pedestal');
end

clear photonimg_dark_tot
in_struc.pedestal=pedestal;

%% Get the vacuum image for subtraction as well:
if ismember(runnum,179:255)
    relevantVacuumRun=178;
    if runnum==179
        run('calc_Xe.m');
        load('pixMap2_temp');
    end
    
    % Run LargeRun for the vacuum run (178) and cache the image.
    [~,~,VacNoUVimg,VacGoodimg,VacIPMdata] = cache_results(@LargeRun,{experiment,relevantVacuumRun,timebins_fs,10000,1}); % caching version %bms, added timebins to LargeRun input (it wasn't working before)
    
    for s=1:length(steps)
        goodipm(s) = mean(IPMandTTdata{s}.goodIPM);
        noUVipm(s)=mean(IPMandTTdata{s}.nouvIPM);
    end
    VacNoUVipm = mean(VacIPMdata{1}.nouvIPM);VacGoodipm=mean(VacIPMdata{1}.goodIPM);
%     modifiedXenon = reshape(modified_all_2D.Xe.*AreaFlat.*lFlat.*BeFlat.*AirFlat.*XeFlat./R2Flat,[388 185 32]);
%     VacNoUVimg = zeros([388 185 32]); VacGoodimg = zeros([388 185 32]); %VacNoUVimg+25*modifiedXenon; VacGoodimg = VacGoodimg+25*modifiedXenon;
elseif ismember(runnum,260:327)
    relevantVacuumRun=259;
    [~,~,VacNoUVimg,VacGoodimg,VacIPMdata] = cache_results(@LargeRun,{experiment,relevantVacuumRun,timebins_fs,10000,1}); % caching version %bms, added timebins to LargeRun input (it wasn't working before)
    for s=1:length(steps)
        goodipm(s)=mean(IPMandTTdata{s}.goodIPM);noUVipm(s)=mean(IPMandTTdata{s}.nouvIPM);
    end
    VacNoUVipm=mean(VacIPMdata{1}.nouvIPM);VacGoodipm=mean(VacIPMdata{1}.goodIPM);
else
    disp('No vacuum found for this run, not subtracting any vacuum!')
    VacNoUVimg=zeros([388 185 32]);
    VacGoodimg=zeros([388 185 32]);
    goodipm=ones([1 length(steps)]);noUVipm=ones([1 length(steps)]);VacNoUVipm=1;VacGoodipm=1;
end
%% No matter what, we're skipping 2d vacuum subtraction
disp('No vacuum found for this run, not subtracting any vacuum!')
    VacNoUVimg=zeros([388 185 32]);
    VacGoodimg=zeros([388 185 32]);
    goodipm=ones([1 length(steps)]);noUVipm=ones([1 length(steps)]);VacNoUVipm=1;VacGoodipm=1;

%% Perform full analysis on each timestep:
% Create empty matrices to house the binned images here so that starting a new step will not wipe them out.
img_good_binned_tot=zeros([388 185 32 length(in_struc.timebins_sec)]); N_good_binned_tot=0; img_good_binned_avg=img_good_binned_tot; IPM_good_binned_tot=0;
img_noUV_binned_tot=zeros([388 185 32 length(in_struc.timebins_sec)]); N_noUV_binned_tot=0; img_noUV_binned_avg=img_noUV_binned_tot; IPM_noUV_binned_tot=0;

for s=steps
    %% Reset variables pertaining to the current timestep:
    disp(['Starting step ' num2str(1+s-steps(1)) ' of ' num2str(length(steps))]);stepstarttime=tic;
%     if Wavelength(s)==0
%         continue
%     end
    in_struc.step=s;
    currentgoodshots=goodshots{s};currentnoUVshots=noUVshots{s};
    in_struc.IPMandTTdata=IPMandTTdata{s};save('tempIPMdata.mat','IPMandTTdata','-v7.3')
    
    % Set current timestep (angle shift) for static points:
    in_struc.current_timestep_sec=timesteps_sec(s);
    if ismember(runnum,258:266) && strcmp(experiment,'i0613')
        longsteptimes=(1e-15)*[5000 -5000 2000 -2000 0 -300 300 10000 -10000];
        in_struc.current_timestep_sec = longsteptimes(runnum - 257);
    elseif ismember(runnum,157:167) && strcmp(experiment,'i0613')
        longsteptimes=(1e-15)*[2000 5000 -2000 -5000 Inf Inf 0  300 -300 600 -600];
        in_struc.current_timestep_sec = longsteptimes(runnum - 156);
    elseif ismember(runnum,66:128) && strcmp(experiment,'b0114')
        longsteptimes=(1e-15)*[-1000 1000 -600 600 -300 300 0 -5000 5000 -2000 2000 2000 -2000 -1000 1000 -600 600 -300 300 0 -5000 5000 5000 -5000 -5000 5000 -2000 2000 -1000 1000 -600 600 -300 300 1000 -600 600 -300 300 0 -5000 5000 -2000 2000 -1000 1000 -600 600 -300 300 0 -5000 5000 -2000 2000 -1000 1000 -600 600 -300 300 0];
        in_struc.current_timestep_sec = longsteptimes(runnum - 65);
    elseif ismember(runnum,[203;206;207;210;212;215]) && strcmp(experiment,'b0114')
        longsteptimes=(1e-15)*2000;
        in_struc.current_timestep_sec = longsteptimes;
    end
    disp(['Timestep identified as ' num2str(1e15*in_struc.current_timestep_sec) ' fs']);
   
    % Ignore extranneous data:
    if in_struc.current_timestep_sec>(max(timebins_sec)+250e-15) || in_struc.current_timestep_sec<(min(timebins_sec)-250e-15);
        disp('Skipping this step because it is outside the desired range.');
        continue
    end
    
    % Reset the shot counters and images from the previous step:
    N_goodshots_counter=0; N_noUVshots_counter=0;
    photonimg_good_unbinned_tot  = zeros([388 185 32]); photonimg_noUV_unbinned_tot = zeros([388 185 32]);
%     img_good_unbinned_tot = zeros([388 185 32]); img_noUV_unbinned_tot = zeros([388 185 32]);
    
    %% Generate new PixMap
    [~,~,q,~,~,pixMap,phivalues,~,~,angleMap,phiMap] = image2rad(experiment,runnum,zeros([388,185,32]),zeros([388,185,32]),D,ones([388,185,32]),Wavelength(s),x0,y0); % Generate new pixMap %bms, edited output

    %% Fetch pressures for this step:
    [pressures_torr] = fetch_pressures(fina,s);
    avg_pressure_torr=mean(pressures_torr);
    
    %% Skip this step if there are no good shots:
    if isempty(currentgoodshots) && isdarkrun==0
        disp('Skipping this step due to missing UVon shots.');
        continue
    elseif isempty(currentnoUVshots)
        disp('WARNING: Missing no-UV data, so all data is treated as such.')
        currentnoUVshots=currentgoodshots;
        in_struc.IPMandTTdata.noUVTTdata_sec=in_struc.IPMandTTdata.goodTTdata_sec;
        in_struc.IPMandTTdata.noUVTTtsMS=in_struc.IPMandTTdata.goodTTtsMS;
        currentgoodshots=currentgoodshots(1:10);
        in_struc.IPMandTTdata.goodTTdata_sec=in_struc.IPMandTTdata.goodTTdata_sec(1:10);
        in_struc.IPMandTTdata.goodTTtsMS=in_struc.IPMandTTdata.goodTTtsMS(1:10);
    end
    
    %% Make mask of pixels with extreme pedestal values:
    pedestalmode=10*mode(round(pedestal(:)/10));%figure(4);subplot(1,2,1);hist(pedestal(:),50);
    if strcmp(experiment,'56012')
        pedestalmask=(pedestal>1000 & pedestal<2000);
    elseif strcmp(experiment,'i0613') || strcmp(experiment,'b0114')
        pedestalmask_new=(pedestal>((pedestalmode*.85)-6) & pedestal<((pedestalmode*1.25)+6));
        pedestalmask = pedestalmask.*pedestalmask_new;
        disp('Updated pedestal mask.');
    end
    flatpedestal=pedestal.*pedestalmask;%figure(4);subplot(1,2,2);hist(flatpedestal(:),50);
    in_struc.pedestalmask=pedestalmask;goodPixels=goodPixels.*pedestalmask;
    disp(['Number of pixels masked by pedestal constraints: ' num2str(100*(2296960-sum(pedestalmask(:)))/2296960) '%']);
    
%     [I_dark,q,~,~] = image2rad(experiment,runnum,photonimg_dark_avg,D,goodPixels,Wavelength(s),x0,y0);
    figure(40); subplot(2,3,1);imagesc(CsPadRearrangeXPP(goodPixels.*photonimg_dark_avg)); title('Image Average: Dark'); axis square; colorbar;
%     subplot(2,3,4);plot(q,I_dark);xlim([plottingqmin plottingqmax]);axis square;title('Radial Average: Dark');
    pause(.5);
        
    %% (2/3) CRUNCH UV-OFF SHOTS %%
        
    N_noUVshots=length(currentnoUVshots); in_struc.N_shots=N_noUVshots; N_loops=ceil(N_noUVshots/shots_per_loop); % Updates the total number of shots.
    noUVstarttime=tic; disp(['Analyzing ' num2str(N_noUVshots) ' UV-off shots']);
    in_struc.isXRayon=1; in_struc.isUVon=0;
    in_struc.count_sums=[];
    
    for j=1:N_loops
        %% Make list of shots to be used for this loop:
        in_struc.minshot=(j-1)*shots_per_loop+1;
        in_struc.maxshot=j*shots_per_loop;
        try
            in_struc.currentshots=currentnoUVshots(in_struc.minshot:in_struc.maxshot);
        catch
            in_struc.maxshot=length(currentnoUVshots);
            in_struc.currentshots=currentnoUVshots(in_struc.minshot:in_struc.maxshot);
        end
        %% Submit this loop's shots to be combined:
        [img_noUV_binned,N_noUV_binned,IPM_noUV_binned,photonimg_noUV_unbinned,~,N_noUV_unbinned,in_struc]=AnalysisFunc(in_struc);
        
        %% End analysis of the current step if it has run out of frames:
        if N_noUV_unbinned == 0
            countsums_noUV_unbinned=cat(1,countsums_noUV_unbinned,in_struc.count_sums);
            break
        end
        %% Add this loop's data to the running total:
        N_noUVshots_counter=N_noUVshots_counter+N_noUV_unbinned;
        photonimg_noUV_unbinned_tot=photonimg_noUV_unbinned_tot+photonimg_noUV_unbinned;
%         img_noUV_unbinned_tot=img_noUV_unbinned_tot+img_noUV_unbinned;
        img_noUV_binned_tot=img_noUV_binned_tot+img_noUV_binned;
        N_noUV_binned_tot=N_noUV_binned_tot+N_noUV_binned;
        IPM_noUV_binned_tot = IPM_noUV_binned_tot+IPM_noUV_binned;
        
        countsums_noUV_unbinned=cat(1,countsums_noUV_unbinned,in_struc.count_sums);
        
        %% Display performance data:
        last_loop=toc(noUVstarttime); hours_elapsed=floor(last_loop/3600); minutes_elapsed=floor((last_loop-3600*hours_elapsed)/60); seconds_elapsed=floor(last_loop-3600*hours_elapsed-60*minutes_elapsed); minutes_remaining=ceil(((N_loops-j)*last_loop/j)/60);
        disp(['Loop ' num2str(j) ' of ' num2str(N_loops) ' finished, ' num2str(hours_elapsed) ':' num2str(minutes_elapsed,'%02.0f') ':' num2str(seconds_elapsed,'%02.0f') ' elapsed, less than ' num2str(minutes_remaining) ' minutes remaining.']);
    end
    
    %% Make hotmask and deadmask
    [deadmask,hotmask,boundarymask] = universalmaskmaker(STDlimit,photonimg_noUV_unbinned_tot/N_noUVshots_counter,pixMap,q,plottingqmin,plottingqmax);
    goodPixels=pedestalmask.*deadmask.*hotmask.*boundarymask;
    if runnum==86
        save('tempgoodPixels','goodPixels');
    elseif runnum==88 || runnum==93
        tempMap=load('tempgoodPixels');
        goodPixels=goodPixels.*tempMap.goodPixels;
    end
%     goodPixels=pedestalmask.*boundarymask;
    tempGoodPixels=pedestalmask.*deadmask.*hotmask;
    disp(['Number of pixels masked by dead/hot mask constraints: ' num2str(100*(2296960-sum(tempGoodPixels(:)))/2296960) '%']);
    
    photonimg_noUV_unbinned_avg=((photonimg_noUV_unbinned_tot/N_noUVshots_counter).*goodPixels)-repmat(VacNoUVimg*noUVipm(s)/VacNoUVipm,[1 1 1 length(N_noUV_unbinned)]);
    img_noUV_binned_avg=(img_noUV_binned_tot./permute(repmat(N_noUV_binned_tot,[1 388 185 32]),[2 3 4 1]))-repmat(VacNoUVimg*noUVipm(s)/VacNoUVipm,[1 1 1 length(N_noUV_binned_tot)]);
    img_noUV_binned_tot_vacuumsub=(img_noUV_binned_tot)-repmat(VacNoUVimg*noUVipm(s)/VacNoUVipm,[1 1 1 length(N_noUV_binned_tot)]).*permute(repmat(N_noUV_binned_tot,[1 388 185 32]),[2 3 4 1]); %bms, added line to vaccuum subtract without averaging
%     IPM_noUV_binned_avg=mean(IPM_noUV_binned_tot,2);
    
    
    %% Plot with Masks
    
    [I_noUV_unbinned,~,q,~,~,~] = image2rad(experiment,runnum,photonimg_noUV_unbinned_avg,zeros([388,185,32]),D,goodPixels,Wavelength(s),x0,y0);
   
    [I_noUV_unbinned_temp,~,q_temp,~,~,~,I_phi_noUV_unbinned,phivalues] = image2rad(experiment,runnum,photonimg_noUV_unbinned_avg,zeros([388,185,32]),D,goodPixels,Wavelength(s),x0,y0);
    
    figure(40); subplot(2,3,2);imagesc(CsPadRearrangeXPP(photonimg_noUV_unbinned_avg.*goodPixels)); title('Image Average: UV off'); axis square; colorbar;
    subplot(2,3,5);plot(q_temp,I_noUV_unbinned_temp);xlim([plottingqmin plottingqmax]);axis square;title('Radial Average: UV off');
    
    pause(.5);
    
    %% (3/3) CRUNCH UV-ON SHOTS %%
    N_goodshots=length(currentgoodshots);in_struc.N_shots=N_goodshots;N_loops=ceil(N_goodshots/shots_per_loop); % Updates the total number of shots.
    goodstarttime=tic;disp(['Analyzing ' num2str(N_goodshots) ' good shots']);
    in_struc.isXRayon=1; in_struc.isUVon=1;
    in_struc.count_sums=[];
    N_good_binned_separate{s}=0;
    
    for j=1:N_loops
        %% Make list of shots to be used for this loop:
        in_struc.minshot=(j-1)*shots_per_loop+1;
        in_struc.maxshot=j*shots_per_loop;
        try
            in_struc.currentshots=currentgoodshots(in_struc.minshot:in_struc.maxshot);
        catch
            in_struc.maxshot=length(currentgoodshots);
            in_struc.currentshots=currentgoodshots(in_struc.minshot:in_struc.maxshot);
        end
        
        %% Submit this loop's shots to be combined:
        [img_good_binned,N_good_binned,IPM_good_binned,photonimg_good_unbinned,~,N_goodshots_unbinned,in_struc]=AnalysisFunc(in_struc);
        if N_goodshots_unbinned == 0 % Ends analysis of the current step if it has run out of frames.
            countsums_good_unbinned=cat(1,countsums_good_unbinned,in_struc.count_sums);
            break
        end
        N_good_binned_separate{s} = N_good_binned_separate{s} + N_good_binned;
        
        %% Add this loop's data to the running total:
        N_goodshots_counter=N_goodshots_counter+N_goodshots_unbinned;
        photonimg_good_unbinned_tot=photonimg_good_unbinned_tot+photonimg_good_unbinned;
%         img_good_unbinned_tot=img_good_unbinned_tot+img_good_unbinned;
        img_good_binned_tot=img_good_binned_tot+img_good_binned;
        N_good_binned_tot=N_good_binned_tot+N_good_binned;
        countsums_binned_tot=countsums_binned_tot+in_struc.binnedcounttot;
        IPM_good_binned_tot=IPM_good_binned_tot+IPM_good_binned;
        
        countsums_binned_avg=countsums_binned_tot./N_good_binned_tot; % Don't divide here anyway
        figure(433);plot(1e15*timebins_sec,countsums_binned_avg);title('Average Total Counts per Timepoint');xlabel('Timebin, fs');ylabel('Counts');
        
        countsums_good_unbinned=cat(1,in_struc.count_sums,countsums_good_unbinned);
        realTimes=cat(1,in_struc.realTimes,realTimes);
        figure(434);scatter(realTimes*1e15,countsums_good_unbinned);title('Detector Intensity Dependence on TT Delay');xlabel('Absolute Time, fs');ylabel('Detector Counts');

        pause(.5);

        %% Display performance data:
        last_loop=toc(goodstarttime); hours_elapsed=floor(last_loop/3600); minutes_elapsed=floor((last_loop-3600*hours_elapsed)/60); seconds_elapsed=floor(last_loop-3600*hours_elapsed-60*minutes_elapsed); minutes_remaining=ceil(((N_loops-j)*last_loop/j)/60);
        disp(['Loop ' num2str(j) ' of ' num2str(N_loops) ' finished, ' num2str(hours_elapsed) ':' num2str(minutes_elapsed,'%02.0f') ':' num2str(seconds_elapsed,'%02.0f') ' elapsed, less than ' num2str(minutes_remaining) ' minutes remaining.']);
    end
    
    photonimg_good_unbinned_avg=((photonimg_good_unbinned_tot/N_goodshots_counter).*goodPixels)-repmat(VacGoodimg*goodipm(1)/VacGoodipm,[1 1 1 length(N_goodshots_unbinned)]);
    img_good_binned_avg=(img_good_binned_tot./permute(repmat(N_good_binned_tot,[1 388 185 32]),[2 3 4 1]))-repmat(VacGoodimg*goodipm(1)/VacGoodipm,[1 1 1 length(N_good_binned_tot)]);
    img_good_binned_tot_vacuumsub=(img_good_binned_tot)-repmat(VacGoodimg*goodipm(s)/VacGoodipm,[1 1 1 length(N_good_binned_tot)]).*permute(repmat(N_good_binned_tot,[1 388 185 32]),[2 3 4 1]); %bms, added line to vacuum subtract without averaging
    [I_good_unbinned,~,q,qerr,goodPixels,qMap,I_phi_good_unbinned,~,pctExcitation,angle_of_q] = image2rad(experiment,runnum,photonimg_good_unbinned_avg,zeros([388,185,32]),D,goodPixels,Wavelength(s),x0,y0); %bms, changed output
    
%     [I_good_unbinned,q,qerr,goodPixels,qMap,I_phi_good_unbinned,~,pctExcitation,angle_of_q] = image2rad(experiment,runnum,photonimg_good_unbinned_avg,D,goodPixels,Wavelength(s),x0,y0);
%     IPM_good_binned_avg=mean(IPM_good_binned_tot,2);
    
    % Plot the last of the data:
    difference_image=goodPixels.*(photonimg_good_unbinned_tot/N_goodshots_counter-photonimg_noUV_unbinned_tot/N_noUVshots_counter);
    figure(40); subplot(2,3,3);imagesc(CsPadRearrangeXPP(difference_image)); title('Image Average Difference'); axis square; colorbar;
    subplot(2,3,6);plot(q,100*(I_good_unbinned-I_noUV_unbinned)./I_noUV_unbinned);xlim([plottingqmin plottingqmax]);axis square;title('Radial Average Difference');legend('%');
    pause(.5);
    
    %% SAVE SINGLE-STEP DATA IF THAT'S ALL THERE IS
    actual_step_position=find(timeorder_steps==s);
    if length(steps)==1
        myfile=['data/',experiment,'/Run',num2str(runnum),'-',num2str(N_goodshots_counter),'frames','-TimeStep',num2str(actual_step_position),'-',my_time_stamp];
        save(myfile,'pctExcitation','x0','y0','Wavelength','D','photonimg_good_unbinned_avg','photonimg_noUV_unbinned_avg','N_goodshots_counter','N_noUVshots_counter','I_good_unbinned','I_noUV_unbinned','q','goodPixels','avg_pressure_torr','IPMandTTdata');
    end
%     save('radial_pixel_histories.mat');
    
    %% ADD CURRENT TIMESTEP TO THE SERIES:
    if s==1
        I_allsteps_good_unbinned=zeros([totalnumsteps length(I_good_unbinned)]);
        I_allsteps_noUV_unbinned=zeros([totalnumsteps length(I_good_unbinned)]);
%         I_allsteps_dark_unbinned=zeros([totalnumsteps length(I_good_unbinned)]);
    end
    
    I_allsteps_good_unbinned(actual_step_position,:)=I_good_unbinned;
    I_allsteps_noUV_unbinned(actual_step_position,:)=I_noUV_unbinned;
%     I_allsteps_dark_unbinned(actual_step_position,:)=I_dark;
    allpressures(actual_step_position)=avg_pressure_torr;
    
    %% Display analysis performance data:
    time_per_step=toc(stepstarttime);
    disp(['Step ' num2str(s) ' took ' num2str(time_per_step/60) ' minutes.']);
    disp(['Good shots had ' num2str(100*mean(countsums_good_unbinned)/mean(countsums_noUV_unbinned)) ' the total detector counts as UV-off']);
end


%% Calculate I for time-binned data:
disp('Generating I for each timepoint.');
if ordered_timepoints_sec==0
    ordered_timepoints_sec=in_struc.current_timestep_sec;
end

% Preallocate matrices:
timebins_sec=in_struc.timebins_sec;
I_good_binned=zeros([length(timebins_sec) length(q)]);
I_good_binned_adu=zeros([length(timebins_sec) length(q)]);
I_noUV_binned=zeros([length(timebins_sec) length(q)]);
I_noUV_binned_adu=zeros([length(timebins_sec) length(q)]);
I_phi_good_binned=zeros([length(timebins_sec) length(phivalues)]);
I_phi_noUV_binned=zeros([length(timebins_sec) length(phivalues)]);


% Calculate I(q) for any bins that actually have data:
for binnumber = 1:length(timebins_sec)
    if N_good_binned_tot(binnumber)~=0
%         [I_good_binned(binnumber,:),q,~,~,~,I_phi_good_binned(binnumber,:),phivalues,~,~] = image2rad(experiment,runnum,img_good_binned_avg(:,:,:,binnumber),D,goodPixels,Wavelength(s),x0,y0);
        [I_good_binned(binnumber,:),I_good_binned_adu(binnumber,:),q,~,~,~,I_phi_good_binned(binnumber,:),phivalues,~,~] = image2rad(experiment,runnum,img_good_binned_avg(:,:,:,binnumber),img_good_binned_tot_vacuumsub(:,:,:,binnumber),D,goodPixels,Wavelength(s),x0,y0); %bms, changed input/output to accommodate new variables
    end
    if N_noUV_binned_tot(binnumber)~=0
%         [I_noUV_binned(binnumber,:),q,~,~,~,I_phi_noUV_binned(binnumber,:),phivalues,~,~,~] = image2rad(experiment,runnum,img_noUV_binned_avg(:,:,:,binnumber),D,goodPixels,Wavelength(s),x0,y0);
        [I_noUV_binned(binnumber,:),I_noUV_binned_adu(binnumber,:),q,~,~,~,I_phi_noUV_binned(binnumber,:),~] = image2rad(experiment,runnum,img_noUV_binned_avg(:,:,:,binnumber),img_noUV_binned_tot_vacuumsub(:,:,:,binnumber),D,goodPixels,Wavelength(s),x0,y0); %bms, changed input/output to accommodate new variables
    end
end

%% GENERATE AND PLOT DIFFERENCE PATTERN
% I_difference_scaled=zeros(size(I_good_binned));
% I_difference=zeros(size(I_good_binned));
% first_existing_point=find(N_good_binned_tot,1,'first');
% for j=1:length(timebins_sec)
%     if N_noUV_binned_tot~=0
%         I_difference_scaled(j,:)=(I_good_binned(j,:)-I_good_binned(first_existing_point,:)).*(q);
%         I_difference(j,:) = I_good_binned(j,:)-I_good_binned(first_existing_point,:);
%     end
% end
%     figure(50);
%     subplot(1,2,2);surf(1e15*timebins_sec(2:end),q(:),I_difference_scaled(2:end,:)');ylim([plottingqmin plottingqmax]);xlim([1e15*min(timebins_sec) 1e15*max(timebins_sec)]);
%     xlabel('Timestep, fs');ylabel('q, inverse Angstrom');title('Percent Change in Scattering');axis square;%shading interp;
%     subplot(1,2,1);surf(1e15*timebins_sec(2:end),q(:),I_difference(2:end,:)');ylim([plottingqmin plottingqmax]);xlim([1e15*min(timebins_sec) 1e15*max(timebins_sec)]);
%     xlabel('Timestep, fs');ylabel('q, inverse Angstrom');title('Absolute Change in Scattering');axis square;%shading interp;

%% Save multi-step data
myfile=['data/',experiment,'/Run',num2str(runnum),'-Allsteps_withIadu-',my_time_stamp];
% save(myfile,'pctExcitation','phivalues','I_phi_good_binned','I_phi_noUV_binned','IPMandTTdata','runnum','countsums_good_unbinned','countsums_noUV_unbinned','countsums_binned_avg','timebins_sec','I_good_binned','N_good_binned_tot','I_noUV_binned','N_noUV_binned_tot','I_allsteps_good_unbinned','I_allsteps_noUV_unbinned','q','qerr','allpressures','ordered_timepoints_sec','N_goodshots_counter','N_noUVshots_counter','N_darkshots_counter','plottingqmax','plottingqmin','N_good_binned_separate');
save(myfile,'pctExcitation','phivalues','I_phi_good_binned','I_phi_noUV_binned','IPMandTTdata','runnum','countsums_good_unbinned','countsums_noUV_unbinned','countsums_binned_avg','timebins_sec','I_good_binned','I_good_binned_adu','N_good_binned_tot','I_noUV_binned','I_noUV_binned_adu','N_noUV_binned_tot','I_allsteps_good_unbinned','I_allsteps_noUV_unbinned','q','qerr','allpressures','ordered_timepoints_sec','N_goodshots_counter','N_noUVshots_counter','N_darkshots_counter','plottingqmax','plottingqmin','N_good_binned_separate','IPM_good_binned_tot','IPM_noUV_binned_tot'); %bms, saved new variables
disp(['Saved as ' myfile]);
disp(['Finished Run ',num2str(runnum)]);
runendtime=toc(runstarttime);
disp(['Run ' num2str(runnum) ' took ' num2str(runendtime/60) ' minutes.']);

end
