function [goodindexes,darkindexes,noUVindexes,IPMandTTdata] = goodshotsonly(hutch,experiment,runnum,step,N_shots)
% This function converts the list of desired shots to only the shots that
% have x-ray intensity as reported by the EVR codes (event codes).
    tic
    fina=GetExpFina(hutch,experiment,runnum);

%% First, the EVR and Time-Tool event codes for each shot (for intentionally dropped shots) are fetched:
    temp=rdXPPdata(fina,'EVRdata','EVRdataTS');
%     temp=rdXPPdata(fina,'EVRdata','EVRdataTS','tt','ttTS');
    
    EVRcode1=zeros(size(temp.scandata(step).EVRdata.rawdata))';
    EVRcode2=zeros(size(temp.scandata(step).EVRdata.rawdata))';
    EVRdataTSms=temp.scandata(step).EVRdata.timeINms;
%     TTdata=temp.scandata(step).tt.FLTPOS_PS;
%     TTtsMS=temp.scandata(step).tt.timeINms;
%     TTampl=temp.scandata(step).tt.AMPL;
    
    for i=1:length(EVRcode1)
        EVRcode1(i)=temp.scandata(step).EVRdata.rawdata(i).eventCode(1);  % Picks the FIRST event code only.
        if strcmp(experiment,'56012')
        else
            EVRcode2(i)=temp.scandata(step).EVRdata.rawdata(i).eventCode(2);
        end
    end
    
    %% Get TT data:
    try
%         [ TTdata,TTampl,TTtsMS,TTpixel ] = fitOpal(fina,step);
        [ TTdata,TTampl,TTtsMS,TTpixel ] = cache_results(@fitOpal,{fina,step});
    catch
        disp('NO TIME-TOOL DATA FOUND!');
        TTdata=zeros(size(EVRcode1));
        TTampl=100*ones(size(TTdata));
        TTtsMS=EVRdataTSms;
        TTpixel=512*ones(size(TTampl));
    end
        
    %% Load CSPAD timestamp and IPM data:
    CSPADtimestr = rdHDFsimple(fina,['/Configure:0000/Run:0000/CalibCycle:' num2str(step-1,'%04d') '/CsPad::ElementV2/XppGon.0:Cspad.0/time']);
    CSPADtsMS = 1000*double(CSPADtimestr.seconds-min(CSPADtimestr.seconds)) + round(1e-5*double(CSPADtimestr.nanoseconds))/10;
    CSPADtsMSmaster=CSPADtsMS;
    try
        [IPMdata,IPMdataTSms,IPMxpos,IPMypos]=GetIPM3Sum(fina,step);
    catch
        disp('Warning: IPM3 not used.');
        [IPMdata,IPMdataTSms,IPMxpos,IPMypos]=GetIPM1Sum(fina,step);
        IPMxpos=ones(size(IPMxpos));
        IPMypos=ones(size(IPMypos));
    end
    
    %% Generate fake numbers for dark runs
    if (strcmp(experiment,'i0613') && runnum==173) || strcmp(experiment,'56012')
        disp('TT and IPM data are ignored for this run.')
        TTdata=zeros(size(TTdata));
        TTampl=100*ones(size(TTampl));
        IPMdata=1e-14*ones(size(IPMdata));
        IPMxpos=ones(size(IPMxpos));
        IPMypos=ones(size(IPMypos));
    end
    
    %% For diagnostics, plot fetched data:
    figure(15);
        subplot(3,2,1);hist(TTdata,50);title('TT jitter for all frames'); xlabel('Picoseconds, uncentered');
        subplot(3,2,2);hist(TTampl,50);title('TT fit amplitudes for all frames');
    figure(16);
        subplot(3,2,1);plot(TTdata);title('TT jitter for all frames'); xlabel('Frame Number');
        subplot(3,2,2);plot(TTampl);title('TT fit amplitudes for all frames');
    pause(.5);
    
    %% Compare the timestamps for CSPAD images, EVR codes, and IPM data and toss any that don't match.
    usableCSPADindexes = find(ismember(CSPADtsMS,EVRdataTSms) & ismember(CSPADtsMS,IPMdataTSms) & ismember(CSPADtsMS,TTtsMS));
    usableEVRindexes   = find(ismember(EVRdataTSms,CSPADtsMS) & ismember(EVRdataTSms,IPMdataTSms) & ismember(EVRdataTSms,TTtsMS));
    usableIPMindexes   = find(ismember(IPMdataTSms,CSPADtsMS) & ismember(IPMdataTSms,EVRdataTSms) & ismember(IPMdataTSms,TTtsMS));
    usableTTindexes    = find(ismember(TTtsMS,CSPADtsMS) & ismember(TTtsMS,IPMdataTSms) & ismember(TTtsMS,EVRdataTSms));
    
    CSPADtsMS=CSPADtsMS(usableCSPADindexes);
    EVRcode1=EVRcode1(usableEVRindexes);
    EVRcode2=EVRcode2(usableEVRindexes);
    EVRdataTSms=EVRdataTSms(usableEVRindexes);
    IPMdata=IPMdata(usableIPMindexes);
    IPMdataTSms=IPMdataTSms(usableIPMindexes);
    IPMxpos=IPMxpos(usableIPMindexes);
    IPMypos=IPMypos(usableIPMindexes);
    TTdata=TTdata(usableTTindexes);
    TTampl=TTampl(usableTTindexes);
    TTtsMS=TTtsMS(usableTTindexes);
    
    if (CSPADtsMS~=EVRdataTSms) | (CSPADtsMS~=IPMdataTSms) | (CSPADtsMS~=TTtsMS)
        disp('Not all timestamps match!');
        save('ERROR-timestamps.mat','usableCSPADindexes','usableEVRindexes','usableIPMindexes','usableTTindexes','CSPADtsMS','EVRdataTSms','IPMdataTSms','TTtsMS');
    end
    
    if isempty(usableCSPADindexes)
        disp('Failed to CSPAD data to other shot data.')
    end
    
    if strcmp(experiment,'i0613')
        %% Assign shots into groups based on EVR codes (note: these should be row, not column vectors)
        %  (Note: For the right CSPAD images to be fetched, it has to come up with an index number out of ALL CSPAD shots, not just the good ones):
        goodindexes=usableCSPADindexes(EVRcode1==140 & EVRcode2~=162)';disp([num2str(length(goodindexes)) ' good shots found.']);
        darkindexes=usableCSPADindexes(EVRcode1==162 | EVRcode1==163 | EVRcode2==162)';disp([num2str(length(darkindexes)) ' dark shots found.']);
        noUVindexes=usableCSPADindexes(EVRcode1==41 | EVRcode1==67 | EVRcode1==141 & EVRcode2~=162)';disp([num2str(length(noUVindexes)) ' noUV shots found.']);
        unidentified_shots=usableCSPADindexes(EVRcode1~=140 & EVRcode1~=162 & EVRcode1~=163 & EVRcode1~=41 & EVRcode1~=67 & EVRcode1~=141)';
        if not(isempty(unidentified_shots))
            disp('Warning: Found unclassified EVR codes.');
        end
        %% Prepare to toss shots by various IPM problems
        goodIPMx=IPMxpos(EVRcode1==140);
        goodIPMy=IPMypos(EVRcode1==140);
        goodIPM=IPMdata(EVRcode1==140);
        nouvIPMx=IPMxpos(EVRcode1==41 | EVRcode1==67 | EVRcode1==141);
        nouvIPMy=IPMypos(EVRcode1==41 | EVRcode1==67 | EVRcode1==141);
        nouvIPM=IPMdata(EVRcode1==41 | EVRcode1==67 | EVRcode1==141);
        goodTTdata=TTdata(EVRcode1==140)-mean(TTdata(EVRcode1==140));
        noUVTTdata=TTdata(EVRcode1==41 | EVRcode1==67 | EVRcode1==141)-mean(TTdata(EVRcode1==41 | EVRcode1==67 | EVRcode1==141));
        goodTTtsMS=TTtsMS(EVRcode1==140);
        noUVTTtsMS=TTtsMS(EVRcode1==41 | EVRcode1==67 | EVRcode1==141);
        goodTTampl=TTampl(EVRcode1==140);
        noUVTTampl=TTampl(EVRcode1==41 | EVRcode1==67 | EVRcode1==141);
        
    elseif strcmp('b0114',experiment)
        goodindexes=usableCSPADindexes(EVRcode1==90 | (EVRcode1==140 & EVRcode2~=162))';disp([num2str(length(goodindexes)) ' good shots found.']);
        noUVindexes=usableCSPADindexes(EVRcode1==91 | EVRcode1==41 & EVRcode2~=162)';disp([num2str(length(noUVindexes)) ' noUV shots found.']);
        darkindexes=usableCSPADindexes(EVRcode1==162 | EVRcode2==162)';disp([num2str(length(darkindexes)) ' dark shots found.']);
        unidentified_shots=usableCSPADindexes(EVRcode1~=90 & EVRcode1~=91 & EVRcode1~=162 & EVRcode1~=40 & EVRcode1~=41)';
        if not(isempty(unidentified_shots))
            disp('Warning: Found unclassified EVR codes.');
        end
        %% Prepare to toss shots by various IPM problems
        goodIPMx=IPMxpos(EVRcode1==90 | EVRcode1==140 & EVRcode2~=162);
        goodIPMy=IPMypos(EVRcode1==90 | EVRcode1==140 & EVRcode2~=162);
        goodIPM=IPMdata(EVRcode1==90 | EVRcode1==140 & EVRcode2~=162);
        nouvIPMx=IPMxpos(EVRcode1==91 | EVRcode1==41 & EVRcode2~=162);
        nouvIPMy=IPMypos(EVRcode1==91 | EVRcode1==41 & EVRcode2~=162);
        nouvIPM=IPMdata(EVRcode1==91 | EVRcode1==41 & EVRcode2~=162);
        TTavg_good=mean(TTdata(EVRcode1==90 | EVRcode1==140 & EVRcode2~=162))
        goodTTdata=TTdata(EVRcode1==90 | EVRcode1==140 & EVRcode2~=162)-TTavg_good;
        TTavg_noUV=mean(TTdata(EVRcode1==91 | EVRcode1==41 & EVRcode2~=162))
        noUVTTdata=TTdata(EVRcode1==91 | EVRcode1==41 & EVRcode2~=162)-TTavg_noUV;
        goodTTtsMS=TTtsMS(EVRcode1==90 | EVRcode1==140 & EVRcode2~=162);
        noUVTTtsMS=TTtsMS(EVRcode1==91 | EVRcode1==41 & EVRcode2~=162);
        goodTTampl=TTampl(EVRcode1==90 | EVRcode1==140 & EVRcode2~=162);
        noUVTTampl=TTampl(EVRcode1==91 | EVRcode1==41 & EVRcode2~=162);
    
    elseif strcmp(experiment,'56012')
        goodindexes=usableCSPADindexes(EVRcode1==90 | EVRcode1==140)';disp([num2str(length(goodindexes)) ' good shots found.']);
        noUVindexes=usableCSPADindexes(EVRcode1==91 | EVRcode1==141)';disp([num2str(length(noUVindexes)) ' noUV shots found.']);
        darkindexes=usableCSPADindexes(EVRcode1==162)';disp([num2str(length(darkindexes)) ' dark shots found.']);
        unidentified_shots=usableCSPADindexes(EVRcode1~=90 & EVRcode1~=91 & EVRcode1~=162 & EVRcode1~=140 & EVRcode1~=141)';
        if not(isempty(unidentified_shots))
            disp('Warning: Found unclassified EVR codes.');
            unique(EVRcode1)
        end
        %% Prepare to toss shots by various IPM problems
        goodIPMx=IPMxpos(EVRcode1==90 | EVRcode1==140);
        goodIPMy=IPMypos(EVRcode1==90 | EVRcode1==140);
        goodIPM=IPMdata(EVRcode1==90 | EVRcode1==140);
        nouvIPMx=IPMxpos(EVRcode1==91 | EVRcode1==141);
        nouvIPMy=IPMypos(EVRcode1==91 | EVRcode1==141);
        nouvIPM=IPMdata(EVRcode1==91 | EVRcode1==141);
        goodTTdata=TTdata(EVRcode1==90 | EVRcode1==140)-mean(TTdata(EVRcode1==90 | EVRcode1==140));
        noUVTTdata=TTdata(EVRcode1==91 | EVRcode1==141)-mean(TTdata(EVRcode1==91 | EVRcode1==141));
        goodTTtsMS=TTtsMS(EVRcode1==90 | EVRcode1==140);
        noUVTTtsMS=TTtsMS(EVRcode1==91 | EVRcode1==141);
        goodTTampl=TTampl(EVRcode1==90 | EVRcode1==140);
        noUVTTampl=TTampl(EVRcode1==91 | EVRcode1==141);
    end

    %% Plot separated TT and IPM data for diagnostics:
    figure(15);
        subplot(3,4,5);hist(goodTTdata,50);title('Good TT data, before filtering');
        subplot(3,4,6);hist(goodTTampl,50);title('Good TT fit amplitude, before filtering');
        subplot(3,4,7);hist(noUVTTdata,50);title('No-UV TT data, before filtering');
        subplot(3,4,8);hist(noUVTTampl,50);title('No-UV TT fit amplitude, before filtering');
    figure(16);
        subplot(3,2,3);plot(1:length(goodTTdata),goodTTdata,1:length(noUVTTdata),noUVTTdata);legend('good TT data','noUV TT data');title('All TT data, before filtering');
        subplot(3,2,4);plot(1:length(goodTTampl),goodTTampl,1:length(noUVTTampl),noUVTTampl);legend('good TT data','noUV TT data');title('All TT fit amplitude, before filtering');
%         subplot(3,4,7);plot(noUVTTdata);title('No-UV TT data, before filtering');
%         subplot(3,4,8);plot(noUVTTampl);title('No-UV TT fit amplitude, before filtering');
    figure(18); %suptitle('X-Ray Centering');
        subplot(2,4,1);scatter(goodIPMx,goodIPMy);title('X-Ray center positions, UV on');%xlim([-.01 .03]);ylim([.01 .05]);axis square;
        subplot(2,4,2);scatter(nouvIPMx,nouvIPMy);title('X-Ray center positions, UV off');%xlim([-.01 .03]);ylim([.01 .05]);axis square;
        subplot(2,4,3);hist(goodIPM,100);title('IPM sum intensity, UV on');axis square;
        subplot(2,4,4);hist(nouvIPM,100);title('IPM sum intensity, UV off');axis square;
        pause(1);
    
    %% Determine undesirable shots (poorly-centered IPM, high or low intensity IPM)
%     averageIPM = mode(round(cat(1,goodIPM,nouvIPM)*10)/10);
    averageIPM = mean(cat(1,goodIPM,nouvIPM));
    centerIPMx_good = mode(round(1000*goodIPMx)/1000); centerIPMx_nouv = mode(round(1000*nouvIPMx)/1000);
    centerIPMy_good = mode(round(1000*goodIPMy)/1000); centerIPMy_nouv = mode(round(1000*nouvIPMy)/1000);
    if strcmp(experiment,'b0114')
        % Generally set at IPM1.25/TT.25 for CHD, IPM1.5/TT.25 for DIB
        uglygoodshots = find(goodIPM<averageIPM/1.25 | goodIPM>averageIPM*1.25 | goodTTampl<.7 | goodTTdata>.25 | goodTTdata<-.25);
        uglynoUVshots = find(nouvIPM<averageIPM/1.25 | nouvIPM>averageIPM*1.25);  
        if ismember(runnum,46:48)
            uglygoodshots = find(goodIPM<averageIPM/1.4 | goodIPM>averageIPM*1.4);
            uglynoUVshots = find(nouvIPM<averageIPM/1.4 | nouvIPM>averageIPM*1.4);
        if runnum==331
            uglygoodshots = find(goodTTampl<.7 | goodTTdata>.25 | goodTTdata<-.25);
            uglynoUVshots = [];
        end
    elseif strcmp(experiment,'i0613')
        uglygoodshots = find(goodIPMy<centerIPMy_good-.0075 | goodIPMy>centerIPMy_good+.0075 | goodIPMx<centerIPMx_good-.001 | goodIPMx>centerIPMx_good+.001 | goodIPM<averageIPM/1.1 | goodIPM>averageIPM*1.1 | goodTTampl<.8 | goodTTdata>.1 | goodTTdata<-.1);
        uglynoUVshots = find(nouvIPMy<centerIPMy_nouv-.0075 | nouvIPMy>centerIPMy_nouv+.0075 | nouvIPMx<centerIPMx_nouv-.001 | nouvIPMx>centerIPMx_nouv+.001 | nouvIPM<averageIPM/1.1 | nouvIPM>averageIPM*1.1);
    elseif strcmp(experiment,'56012')
        uglygoodshots = [];
        uglynoUVshots = [];
    end
    
    num_nouvshots_ugly=length(uglynoUVshots);
    num_goodshots_ugly=length(uglygoodshots);
    num_tot_goodshots=length(goodindexes);
    num_tot_nouvshots=length(noUVindexes);
    
    %% Chuck all data for poorly centered shots
    
    if length(goodIPMx)~=length(goodIPMy) || length(goodIPMx)~=length(goodindexes)
        disp('Pointing index mismatch for good shots.');
    end
    if length(nouvIPMx)~=length(nouvIPMy) || length(nouvIPMx)~=length(noUVindexes)
        disp('Pointing index mismatch for noUV shots.');
    end

    goodIPM(uglygoodshots) = [];
    goodIPMx(uglygoodshots) = [];
    goodIPMy(uglygoodshots) = [];
    nouvIPM(uglynoUVshots) = [];
    nouvIPMx(uglynoUVshots) = [];
    nouvIPMy(uglynoUVshots) = [];
    goodindexes(uglygoodshots) = [];
    noUVindexes(uglynoUVshots) = [];
    goodTTdata(uglygoodshots) = [];
    goodTTampl(uglygoodshots) = [];
    goodTTtsMS(uglygoodshots) = [];
    noUVTTdata(uglynoUVshots) = [];
    noUVTTampl(uglynoUVshots) = [];
    noUVTTtsMS(uglynoUVshots) = [];
    
    if length(goodTTtsMS)==length(CSPADtsMSmaster(goodindexes))
        if goodTTtsMS~=CSPADtsMSmaster(goodindexes)
            disp('Epic fail.')
        end
    end    

    disp([num2str(num_goodshots_ugly) ' UV-on  shots discarded based on IPM or TT. (' num2str(100*num_goodshots_ugly/num_tot_goodshots) '%)']);
    disp([num2str(num_nouvshots_ugly) ' UV-off shots discarded based on IPM or TT. (' num2str(100*num_nouvshots_ugly/num_tot_nouvshots) '%)']);
    
    %% Plot data again for sorted, acceptable data:
    figure(15);
        subplot(3,4,9);hist(goodTTdata,50);title('Good TT data, after filtering');subplot(3,4,10);hist(goodTTampl,50);title('Good TT fit amplitude, after filtering');
        try subplot(3,4,11);hist(noUVTTdata,50);title('No-UV TT data, after filtering');subplot(3,4,12);hist(noUVTTampl,50);title('No-UV TT fit amplitude, before filtering');end;
    figure(16);
        subplot(3,2,5);plot(1:length(goodTTdata),goodTTdata,1:length(noUVTTdata),noUVTTdata);legend('good TT data','noUV TT data');title('All TT data, after filtering');
        subplot(3,2,6);plot(1:length(goodTTampl),goodTTampl,1:length(noUVTTampl),noUVTTampl);legend('good TT data','noUV TT data');title('All TT fit amplitude, after filtering');
%         subplot(3,4,9);plot(goodTTdata);title('Good TT data, after filtering');subplot(3,4,10);plot(goodTTampl);title('Good TT fit amplitude, after filtering');
%         try subplot(3,4,11);plot(noUVTTdata);title('No-UV TT data, after filtering');subplot(3,4,12);plot(noUVTTampl,50);title('No-UV TT fit amplitude, after filtering');end;
    figure(18);
        subplot(2,4,5);scatter(goodIPMx,goodIPMy);title('X-Ray center positions, UV on');%xlim([-.01 .03]);ylim([.01 .05]);axis square;
        try subplot(2,4,6);scatter(nouvIPMx,nouvIPMy);title('X-Ray center positions, UV off');end;%xlim([-.01 .03]);ylim([.01 .05]);axis square;
        subplot(2,4,7);hist(goodIPM,100);title('IPM sum intensity, UV on');axis square;
        try subplot(2,4,8);hist(nouvIPM,100);title('IPM sum intensity, UV off');axis square;end;
        pause(.5);
    
    %% In case you are not analyzing all shots, it only takes the desired number.
    try 
        goodindexes=goodindexes(1:N_shots);
        noUVindexes=noUVindexes(1:N_shots);
        darkindexes=darkindexes(1:N_shots);
    end
    
    %% Rename IPM and TT data for export
    IPMandTTdata.IPMdata=IPMdata; IPMandTTdata.IPMdataTSms=IPMdataTSms; IPMandTTdata.TTpixel=TTpixel;IPMandTTdata.goodIPM=goodIPM;IPMandTTdata.nouvIPM=nouvIPM;
    IPMandTTdata.goodTTdata_sec=(1e-12)*(goodTTdata); IPMandTTdata.goodTTtsMS=goodTTtsMS; IPMandTTdata.TTampl=goodTTampl;
    IPMandTTdata.noUVTTdata_sec=(1e-12)*(noUVTTdata); IPMandTTdata.noUVTTtsMS=noUVTTtsMS;
    
    if (100*num_goodshots_ugly/num_tot_goodshots)>80 || 100*num_nouvshots_ugly/num_tot_nouvshots>80
        goodindexes=[];noUVindexes=[];darkindexes=[];
        disp(['Ignoring data from step ' num2str(step) ' due to many ugly shots.']);
    end
    
%     save('test_goodshot_data.mat')
    toc
end
