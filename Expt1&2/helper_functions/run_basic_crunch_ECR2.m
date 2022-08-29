function [dat,info]=run_basic_crunch_ECR2(dir1)

% for force and position etc. variability comparison

% window size of 225*2+1
includeft = 0; % include force transducer; makes files far bigger.
dbstop if error
info.base_dir = dir1;
info.FFMAG=15;
tic;
d=dir([info.base_dir,'\*.mat']);

wsize=225;
winpad = nan(wsize,1);
count=0;
for k=1:length(d),
    if d(k).isdir==0,
        count=count+1;
        filename_list{count}=d(k).name;
        filename_list_subjectid{count}=d(k).name(1:3);
        filename_list_blockid{count}=d(k).name(4);
    end
end
info.sublist = unique(filename_list_subjectid);
info.blockletter = unique(filename_list_blockid)';
info.blocksize = nan(size(info.blockletter));

dt=dir([info.base_dir,'\*.tgt']);
% if length(dt) > 0,
if ~isempty(dt),
    tgt_filename_list = {dt(:).name};
    for k=1:length(dt),
        tgt_basename{k} = dt(k).name(1:end-5);
        info.blockletter{k} = dt(k).name(end-4);
    end
    tgt_basename = unique(tgt_basename);
    info.blockletter = sort(unique(info.blockletter));
    
    if length(tgt_basename)>1, disp('Error: More than one Tgt file basename detected...'); end
    tgt_basename = tgt_basename{1};
    
    tgt_all = [];
    for k2=1:length(info.blockletter),
        tgt = load([info.base_dir,'\',tgt_basename,info.blockletter{k2},'.tgt']);
        tgt_all = [tgt_all; tgt];
        info.blocksize(k2) = size(tgt,1);
    end
else
    for k1=1:length(info.sublist),
        for k2=1:length(info.blockletter),
            load([info.base_dir,'\',info.sublist{k1},info.blockletter{k2},'.mat'])
            if exist('tgt'),
                if min(size(tgt))==1, tgt = reshape(tgt,8,[])'; end
                tgt_dat{k1,k2}=tgt;
                tgt_dat_length(k1,k2) = size(tgt,1);
            end
        end
    end
    info.blocksize = mode(tgt_dat_length,1)';
    tgt_all = vertcat(tgt_dat{1,:});
    for k2=1:length(info.blockletter),
        for k1a=1:length(info.sublist), for k1b=1:length(info.sublist),
                if length(tgt_dat{k1a,k2}(:)) == length(tgt_dat{k1b,k2}(:)),
                    tgt_dat_equal(k1a,k1b,k2) = all(tgt_dat{k1a,k2}(:) == tgt_dat{k1b,k2}(:));
                else tgt_dat_equal(k1a,k1b,k2) = 0;
                end
                
            end; end
    end
    tgt_dat_equal_summary = all(reshape(tgt_dat_equal,k1a*k1b,k2))';
    if all(tgt_dat_equal_summary), disp('all target files identical (good)!')
    else disp(['Target files for sets ',num2str(find(tgt_dat_equal_summary==0)), ' are NOT identical!! - please check on this'])
    end
    info.tgt_dat = tgt_dat;
    info.tgt_dat_equal = tgt_dat_equal;
    info.tgt_dat_equal_summary = tgt_dat_equal_summary;
    %info.blocksize = tgt_dat_length;
end
info.tgt_all = tgt_all;
%keyboard

v_th = 0.05;
v_tht = 40;
t_startoffset=20;
t_stopoffset=20;
TARGET_DIST=0.1;
n1=10;

display(info)

d=dir([info.base_dir,'\*.mat']);
count=0;
for k=1:length(d),
    if d(k).isdir==0,
        count=count+1;
        %filename_list{count}=d(k).name;
        filename_list_subjectid{count}=d(k).name(1:3);
        %filename_list_blockid{count}=d(k).name(4);
    end
end
info.sublist = unique(filename_list_subjectid);

nm = sum(info.blocksize);
dat.fyr = cell(nm,1);
dat.vxr = cell(nm,1);
dat.vyr = cell(nm,1);
dat.pxr = cell(nm,1);
dat.pyr = cell(nm,1);
dat.ms = cell(nm,1);
if includeft % note: fz contains no data (all zeros)
    dat.fy_ft = cell(nm,1);
    dat.fx_ft = cell(nm,1);
end

for kb = 1:length(info.blocksize), disp(['Block ',info.blockletter(kb)]); % block counter
    for ks = 1:length(info.sublist), disp([info.sublist{ks},' ',num2str(toc)]) % subject counter
        filename = [info.base_dir,'\',info.sublist{ks},info.blockletter{kb},'.mat'];
        q = load(filename);
        if 1, %~isfield(q,'tgt'),
            cs = [0;cumsum(info.blocksize)];
            q.tgt = tgt_all(cs(kb)+1:cs(kb+1),:);
        end
        if min(size(q.tgt))==1, q.tgt = reshape(q.tgt,8,[])'; end;
        if size(q.tgt,1)~=info.blocksize(kb),
            disp('unexpected tgt fle length');
            keyboard
        end
        if max(q.tgt(:,1))<=32, NUM_DIRS=32; else NUM_DIRS=360; end
        phi = -q.tgt(:,1)*2*pi/NUM_DIRS;
        x0_tgt = [0; cumsum(TARGET_DIST*cos(-phi))];
        y0_tgt = [0; cumsum(TARGET_DIST*sin(-phi))];
        % kn is overall movement number)
        for kn1=1:info.blocksize(kb),
            kn=kn1+sum(info.blocksize(1:kb-1));
            if length(q.rawdata)>=kn1,
                if ~isempty(q.rawdata{kn1})
                    x = q.rawdata{kn1}(n1:end,1); y = q.rawdata{kn1}(n1:end,2);
                    vx = q.rawdata{kn1}(n1:end,3); vy = q.rawdata{kn1}(n1:end,4);
                    ax = q.rawdata{kn1}(n1:end,5); ay = q.rawdata{kn1}(n1:end,6);
                    fx = q.rawdata{kn1}(n1:end,9); fy = q.rawdata{kn1}(n1:end,10);
                    timing_data = q.rawdata{kn1}(n1:end,13);
                    
                    %Andrew added these 4 lines to compute the movement
                    %times according to the program.  That way, you can be
                    %sure of the times used by the robot to compute
                    %rewards.
                    mov_stage = q.rawdata{kn1}(n1:end,14);
                    irew = find(q.rawdata{kn1}(:,14)==14);
                    if length(irew)~=1,
                        display('Error: more than one iteration with movement_stage = 14')
                        irew = irew(1);
                        keyboard
                    end
                    dat.rew_code(kn,ks) = q.rawdata{kn1}(irew+1,15);
                    dat.online_AC(kn,ks) = q.rawdata{kn1}(irew+1,16);
                    
                    t_last = find(mov_stage==4,1,'last');
                    t_first = find(mov_stage==4,1,'first')-1;
                    if ~isempty(t_first+t_last)
                        dat.mtime_real(kn,ks) = (t_last-t_first)*0.005;
                    else
                        dat.mtime_real(kn,ks) = nan;
                    end
                    % buffer fx and fy
                    % firstnonzero=find(fx(1:end-10)~=0,1);
                    % fx(1:firstnonzero)=nanmean(fx((firstnonzero+1):(firstnonzero+10)));
                    % firstnonzero=find(fy(1:end-10)~=0,1);
                    % fy(1:firstnonzero)=nanmean(fy((firstnonzero+1):(firstnonzero+10)));
                    fx(mov_stage<2 | mov_stage==7)=nan;
                    fy(mov_stage<2 | mov_stage==7)=nan;
                    
                    fx_ft = q.rawdata{kn1}(n1:end,7); fy_ft = q.rawdata{kn1}(n1:end,8);
                    %[xr,yr,vxr,vyr] = rotate_data(q);
                    
                    % if exist('xr','var'), xprevious = xr; else xprevious = []; end
                    
                    if isfield(q,'tgt'),
                        [xr,yr] = rotate2(phi(kn1), x-x0_tgt(kn1), y-y0_tgt(kn1));
                        [vxr,vyr] = rotate2(phi(kn1), vx, vy);
                        [axr,ayr] = rotate2(phi(kn1), ax, ay);
                        [fxr,fyr] = rotate2(phi(kn1), fx, fy);
                        [fxr_ft,fyr_ft] = rotate2(phi(kn1), fx_ft, fy_ft);
                    else %assume 90/270
                        if rem(kn1,2)==1,
                            xr=-y; yr=x; vxr=-vy; vyr=vx; axr=-ay; ayr=ax; fxr=-fy; fyr=fx; fxr_ft=-fy_ft; fyr_ft=fx_ft; %odd mvts (down)
                        else
                            xr=y+0.1; yr=-x; vxr=vy; vyr=-vx; axr=ay; ayr=-ax; fxr=fy; fyr=-fx; fxr_ft=fy_ft; fyr_ft=-fx_ft; %even mvts (up)
                        end
                    end
                    if phi>(pi/2+.0001) & phi < (1.5*pi-0.001),
                        keyboard;
                    end
                    
                    % t40 is last time at which movement passes 4 cm
                    %t40 = max(find([[xr(1:end-1)<=0.04]&[xr(2:end)>0.04]]));
                    t40 = find((xr(1:end-1)<=0.04) & (xr(2:end)>0.04),1,'last');
                    if t40 >200,
                        vxr(1:t40-200)=nanmean(vxr((t40-200):(t40-190)));
                        vyr(1:t40-200)=nanmean(vyr((t40-200):(t40-190)));
                    end
                    
                    vt = sqrt(vxr.^2+vyr.^2);
                    
                    % Determine whether they jumped the gun or made a pre-movement (before the system even displayed the target and during which time the EC or FF is not on)
                    % Jumped the gun is computed as tanvel>=K_VEL_THRESHOLD_START during show_tgt_nogo, which is movement stage 2.
                    % Also, the requirement to start the next movement is (tanvel >= K_VEL_THRESHOLD_END) || (origin_dist > K_TARGET_START_RADIUS)
                    %    for (staystilltime >= K_STARTMVT_STAYSTILLPERIOD) during move_to_correct_starting_pos, which is movement stage 2.
                    %    These are the cases where they have moved without the FF.
                    ms = q.rawdata{kn1}(n1:end,14); % movement stage
                    
                    % Flag the jump-the-guns
                    last_ms2 = find(ms==2,1,'last');
                    if ~isempty(last_ms2) && length(vt)>last_ms2 && (vt(last_ms2+1)>v_th),
                        dat.jtg(kn,ks) = 1;
                        shg1 = 1;
                    else
                        dat.jtg(kn,ks) = 0;
                        shg1 = 0;
                    end
                    
                    % Flag the pre-movements
                    first_ms1 = find(ms==1,1,'first');
                    last_ms1 = find(ms==1,1,'last');
                    if isempty(first_ms1), first_ms1 = 1; end
                    if isempty(last_ms1), last_ms1 = 1; end
                    if ( (max(xr(first_ms1:last_ms1)-mean(xr(1:10)))>0.01 && max(xr(first_ms1:last_ms1))>0.02) || max(xr(first_ms1:last_ms1))>0.03)  && kn1>1,
                        dat.pm(kn,ks) = 1;
                        shg3 = 1;
                    else
                        shg3 = 0;
                        dat.pm(kn,ks) = 0;
                    end
                    % Check for things not flagged as pre-movements that might qualify
                    if max(xr(first_ms1:last_ms1))>0.02 && shg3~=1 && kn1>1, shg4 = 1; else shg4 = 0; end
                    % Check for lateral movements that might be problematic
                    %if max(abs(yr(1:last_ms2)))>0.02 && kn1>1 && shg3==0, shg5=1; else shg5=0; end
                    
                    if 0 %shg1
                        %display([shg1 shg3 shg4])
                        min1 = find(ms==1,1,'first');
                        max1 = find(ms==1,1,'last');
                        min2 = find(ms==2,1,'first');
                        max2 = find(ms==2,1,'last');
                        if ~isempty(find(ms(1:max2)>2)), display('non-monotone ms'), keyboard, end
                        figure(1), hold off,
                        fill([min1 max1 max1 min1],[0 0 1 1],[1 1 1]*0.9), hold on
                        fill([min2 max2 max2 min2],[0 0 1 1],[1 1 1]*0.7)
                        plot(vt,'r')
                        plot(xr,'b')
                        ylim([-0.03 0.2])
                        display(kn)
                        figure(2), hold off
                        plot(xprevious)
                        figure(3), hold off
                        plot(yr)
                    end
                    
                    % qra is number of samples where person is moving faster
                    % than threshold
                    qra = cumsum(vt>v_th); % v_th is threshold for movement (.05)
                    %ra=min(find([qra(v_tht+1:end)-qra(1:end-v_tht)==v_tht]));
                    ra=find(qra(v_tht+1:end)-qra(1:end-v_tht)==v_tht,1,'first');
                    % ra is last time the person is under threshold
                    
                    % compute ITI:
                    [~,im] = max(abs(vxr));
                    time_before = sum(q.rawdata{kn1}(1:(n1+im-1),13));
                    if kn1==1,
                        dat.iti(kn,ks) = nan;
                        %dat.time_after(kn,ks) = time_after;
                    else
                        %dat.time_before(kn,ks) = time_before;
                        dat.iti(kn,ks) = time_after + time_before;
                        %dat.time_after(kn,ks) = time_after;
                    end
                    time_after = sum(q.rawdata{kn1}((n1+im):end,13));
                    
                    % if person is begins too early or too late then it's not good
                    if isempty(ra) || ra>7000, dat.good(kn,ks) = 0; message='start late';
                    elseif ra<20, dat.good(kn,ks) = 0; message='start early';
                    elseif isempty(t40), dat.good(kn,ks) = 0; message='miss 4cm';
                    elseif dat.iti(kn,ks)>(info.tgt_all(kn,end)+10), dat.good(kn,ks) = 0; message='waited >10s before starting';
                    else dat.good(kn,ks) = 1;
                    end
                    
                    % if movement is too fast or too slow (velthresh)
                    % or if it takes too long (mtime_real) (Andrew added,
                    % Jan 2016).
                    vmaxthresh=0.8;%.5;
                    vminthresh=0.2;%.23;
                    mtmaxthresh=2;%
                    dat.timing_data{kn,ks}=timing_data;
                    %dat.movement_stage{kn,ks} = q.rawdata{kn1}(n1:end,14);
                    dat.vmax(kn,ks)=max(abs(vxr));
                    if max(abs(vxr))>vmaxthresh || max(abs(vxr))<vminthresh || dat.mtime_real(kn,ks)>mtmaxthresh
                        dat.good(kn,ks)=0;
                        if max(abs(vxr))>vmaxthresh, message='too fast';
                        elseif max(abs(vxr))<vminthresh, message='too slow';
                        else message='mvmt took too long';
                        end
                    end
                    
                    %                    figure(99); subplot(211); plot([vxr,vyr]); title(num2str([kn,kn1,ks,ra,t40])); subplot(212); plot([xr,yr]); shg; pause;
                    if  dat.good(kn,ks),
                        
                        % t0 where the actual movement starts
                        t0 = max(1,ra-t_startoffset);
                        % tvtmax = min(find(vt==max(vt))); % time of first max vel
                        % tvxmax = min(find(abs(vxr)==max(abs(vxr))));
                        % tvymax = min(find(abs(vyr)==max(abs(vyr))));
                        
                        tvtmax = find(vt==max(vt),1,'first'); % time of first max vel
                        tvxmax = find(abs(vxr)==max(abs(vxr)),1,'first');
                        tvymax = find(abs(vyr)==max(abs(vyr)),1,'first');
                        
                        tf = min(length(vxr), t0+2*(tvtmax-t0));
                        % end of movement, either the total length or twice the
                        % time from t0 to tvtmax after t0
                        
                        dat.mtime(kn,ks) = tf-t0; %tstop-ra; %tf-t0;
                        dat.pathlen(kn,ks) = sum(vt)*5;
                        dat.mdist(kn,ks) = sqrt((xr(tf)-xr(t0))^2+(yr(tf)-yr(t0))^2);
                        dat.pathlen_rel(kn,ks) = dat.pathlen(kn,ks)/dat.mdist(kn,ks);
                        dat.vtmax(kn,ks) = max(vt);
                        dat.vxmax(kn,ks) = vxr(tvxmax);
                        dat.vymax(kn,ks) = vyr(tvymax);
                        dat.ra(kn,ks) = ra;
                        dat.tvtmax(kn,ks) = tvtmax - ra;
                        dat.tvxmax(kn,ks) = tvxmax - ra;
                        dat.n(kn,ks) = length(vxr);
                        
                        % Andrew's changes: pass on all the data, with nans
                        % where it is missing.
                        
                        % First, pad all vectors with 225 nans before and
                        % after.  Then take the 225 values around maxvel;
                        % this way, if there are not 225 values, it will
                        % grab from the nan padding I added.
                        fyr_pad = [winpad; fyr; winpad];
                        vxr_pad = [winpad; vxr; winpad];
                        vyr_pad = [winpad; vyr; winpad];
                        xr_pad = [winpad; xr; winpad];
                        yr_pad = [winpad; yr; winpad];
                        ms_pad = [winpad; mov_stage; winpad];
                        
                        win = (-wsize:wsize)+tvxmax+wsize;
                        dat.vxr{kn}(:,ks) = vxr_pad(win);
                        dat.vyr{kn}(:,ks) = vyr_pad(win);
                        dat.pxr{kn}(:,ks) = xr_pad(win);
                        dat.pyr{kn}(:,ks) = yr_pad(win);
                        dat.ms{kn}(:,ks)  = ms_pad(win);
                        if q.tgt(kn1,7),
                            dat.fyr{kn}(:,ks) = fyr_pad(win);
                        end
                        if includeft
                            fy_ft_pad = [winpad; fy_ft; winpad];
                            fx_ft_pad = [winpad; fx_ft; winpad];
                            dat.fy_ft{kn}(:,ks) = fy_ft_pad(win);
                            dat.fx_ft{kn}(:,ks) = fx_ft_pad(win);
                        end
                        
                        
                        % if q.tgt(kn1,7),               %if FC trial then save force data
                        %     if tvxmax > wsize, % time of max x velocity is more than windowsize
                        %         if tvxmax+wsize > length(fyr)
                        %             dat.fyr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             dat.vxr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             dat.vyr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             dat.pxr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             dat.pyr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             dat.ms{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %
                        %         else
                        %             dat.fyr{kn}(:,ks) = fyr(tvxmax + [-wsize:wsize]);
                        %             dat.vxr{kn}(:,ks) = vxr(tvxmax + [-wsize:wsize]);
                        %             dat.vyr{kn}(:,ks) = vyr(tvxmax + [-wsize:wsize]);
                        %
                        %             %                            if mod(kn,2)==1
                        %             dat.pxr{kn}(:,ks) = xr(tvxmax + [-wsize:wsize]);
                        %             dat.pyr{kn}(:,ks) = yr(tvxmax + [-wsize:wsize]);
                        %             dat.ms{kn}(:,ks) = mov_stage(tvxmax + [-wsize:wsize]);
                        %             %                            else
                        %             %                                dat.pxr{kn}(:,ks) = -xr(tvxmax + [-wsize:wsize])+.1;
                        %             %                            end
                        %             if includeft
                        %                 dat.fy_ft{kn}(:,ks) = fyr_ft(tvxmax + [-wsize:wsize]);
                        %                 dat.fx_ft{kn}(:,ks) = fxr_ft(tvxmax + [-wsize:wsize]);
                        %             end
                        %         end
                        %     else
                        %         if tvxmax+wsize > length(fyr)
                        %             dat.fyr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             dat.vxr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             dat.pxr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             dat.vyr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             dat.pyr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             dat.ms{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             if includeft
                        %                 dat.fy_ft{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %                 dat.fx_ft{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             end
                        %         else
                        %             % pad vxr and fyr in front
                        %             padsize=wsize-tvxmax+1;
                        %             padvalfyr=nanmean(fyr(find(fyr~=0,10)));
                        %             padvalvxr=nanmean(vxr(find(vxr~=0,10)));
                        %             padvalvyr=nanmean(vyr(find(vyr~=0,10)));
                        %             %                                if mod(kn,2)==1
                        %             padvalxr=nanmean(xr(find(xr~=0,10)));
                        %             padvalyr=nanmean(yr(find(yr~=0,10)));
                        %             %                                else
                        %             %                                    padvalxr=-nanmean(xr(1:10))+.1;
                        %             %                                end
                        %             dat.fyr{kn}(1:padsize,ks)=padvalfyr*ones(padsize,1);
                        %             dat.vxr{kn}(1:padsize,ks)=padvalvxr*ones(padsize,1);
                        %             dat.pxr{kn}(1:padsize,ks)=padvalxr*ones(padsize,1);
                        %             dat.vyr{kn}(1:padsize,ks)=padvalvyr*ones(padsize,1);
                        %             dat.pyr{kn}(1:padsize,ks)=padvalyr*ones(padsize,1);
                        %             dat.ms{kn}(1:padsize,ks) = nan(padsize,1,class(fx));
                        %
                        %             dat.fyr{kn}((padsize+1):(wsize*2+1),ks) = fyr(1:(tvxmax + wsize));
                        %             dat.vxr{kn}((padsize+1):(wsize*2+1),ks) = vxr(1:(tvxmax + wsize));
                        %             if includeft
                        %                 padvalfy_ft=nanmean(fyr_ft(find(fyr_ft~=0,10)));
                        %                 padvalfx_ft=nanmean(fxr_ft(find(fxr_ft~=0,10)));
                        %
                        %                 dat.fy_ft{kn}(1:padsize,ks)=padvalfy_ft*ones(padsize,1);
                        %                 dat.fx_ft{kn}(1:padsize,ks)=padvalfx_ft*ones(padsize,1);
                        %
                        %                 dat.fy_ft{kn}((padsize+1):(wsize*2+1),ks) = fyr_ft(1:(tvxmax + wsize));
                        %                 dat.fx_ft{kn}((padsize+1):(wsize*2+1),ks) = fxr_ft(1:(tvxmax + wsize));
                        %             end
                        %             %                                if mod(kn,2)==1
                        %             dat.vyr{kn}((padsize+1):(wsize*2+1),ks) = vyr(1:(tvxmax + wsize));
                        %             dat.pxr{kn}(padsize+1:wsize*2+1,ks) = xr(1:tvxmax + wsize);
                        %             dat.pyr{kn}(padsize+1:wsize*2+1,ks) = yr(1:tvxmax + wsize);
                        %             dat.ms{kn}(padsize+1:wsize*2+1,ks) = mov_stage(1:tvxmax + wsize);
                        %             %                                else
                        %             %                                    dat.pxr{kn}(padsize+1:wsize*2+1,ks) = -xr(1:tvxmax + wsize)+.1;
                        %             %                                end
                        %         end
                        %     end
                        % else
                        %     % save velocity and position data if non error clamp
                        %     if tvxmax > wsize, % time of max x velocity is more than windowsize
                        %         if tvxmax+wsize > length(vxr)
                        %             dat.vxr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             dat.vyr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             dat.pxr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             dat.pyr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             dat.ms{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             if includeft
                        %                 dat.fy_ft{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %                 dat.fx_ft{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             end
                        %         else
                        %             dat.vxr{kn}(:,ks) = vxr(tvxmax + [-wsize:wsize]);
                        %             dat.vyr{kn}(:,ks) = vyr(tvxmax + [-wsize:wsize]);
                        %
                        %             %                            if mod(kn,2)==1
                        %             dat.pxr{kn}(:,ks) = xr(tvxmax + [-wsize:wsize]);
                        %             dat.pyr{kn}(:,ks) = yr(tvxmax + [-wsize:wsize]);
                        %             dat.ms{kn}(:,ks) = mov_stage(tvxmax + [-wsize:wsize]);
                        %             %                            else
                        %             %                                dat.pxr{kn}(:,ks) = -xr(tvxmax + [-wsize:wsize])+.1;
                        %             %                            end
                        %             if includeft
                        %                 dat.fy_ft{kn}(:,ks)=fyr_ft(tvxmax + [-wsize:wsize]);
                        %                 dat.fx_ft{kn}(:,ks)=fxr_ft(tvxmax + [-wsize:wsize]);
                        %             end
                        %         end
                        %     else
                        %         if tvxmax+wsize > length(vxr)
                        %             dat.vxr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             dat.pxr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             dat.vyr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             dat.pyr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             dat.ms{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             if includeft
                        %                 dat.fy_ft{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %                 dat.fx_ft{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %             end
                        %         else
                        %             % pad vxr and fyr in front
                        %             padsize=wsize-tvxmax+1;
                        %             padvalvxr=nanmean(vxr(find(vxr~=0,10)));
                        %             padvalvyr=nanmean(vyr(find(vyr~=0,10)));
                        %             %                                if mod(kn,2)==1
                        %             padvalxr=nanmean(xr(find(xr~=0,10)));
                        %             padvalyr=nanmean(yr(find(yr~=0,10)));
                        %             %                                else
                        %             %                                    padvalxr=-nanmean(xr(1:10))+.1;
                        %             %                                end
                        %             dat.vxr{kn}(1:padsize,ks)=padvalvxr*ones(padsize,1);
                        %             dat.pxr{kn}(1:padsize,ks)=padvalxr*ones(padsize,1);
                        %             dat.vyr{kn}(1:padsize,ks)=padvalvyr*ones(padsize,1);
                        %             dat.pyr{kn}(1:padsize,ks)=padvalyr*ones(padsize,1);
                        %             dat.ms{kn}(1:padsize,ks) = nan(padsize,1,class(fx));
                        %
                        %             dat.vxr{kn}((padsize+1):(wsize*2+1),ks) = vxr(1:(tvxmax + wsize));
                        %             dat.vyr{kn}((padsize+1):(wsize*2+1),ks) = vyr(1:(tvxmax + wsize));
                        %
                        %             %                                if mod(kn,2)==1
                        %             dat.pxr{kn}(padsize+1:wsize*2+1,ks) = xr(1:tvxmax + wsize);
                        %             dat.pyr{kn}(padsize+1:wsize*2+1,ks) = yr(1:tvxmax + wsize);
                        %             dat.ms{kn}(padsize+1:wsize*2+1,ks) = mov_stage(1:tvxmax + wsize);
                        %             %                                else
                        %             %                                    dat.pxr{kn}(padsize+1:wsize*2+1,ks) = -xr(1:tvxmax + wsize)+.1;
                        %             %                                end
                        %             if includeft
                        %                 padvalfy_ft=nanmean(fyr_ft(find(fyr_ft~=0,10)));
                        %                 padvalfx_ft=nanmean(fxr_ft(find(fxr_ft~=0,10)));
                        %
                        %                 dat.fy_ft{kn}(1:padsize,ks)=padvalfy_ft*ones(padsize,1);
                        %                 dat.fx_ft{kn}(1:padsize,ks)=padvalfx_ft*ones(padsize,1);
                        %
                        %                 dat.fy_ft{kn}((padsize+1):(wsize*2+1),ks) = fyr_ft(1:(tvxmax + wsize));
                        %                 dat.fx_ft{kn}((padsize+1):(wsize*2+1),ks) = fxr_ft(1:(tvxmax + wsize));
                        %             end
                        %         end
                        %     end
                        % end
                    else disp(['aborted (badtrial) #', num2str(kn),': ',message])
                        nanfill = nan(2*wsize+1,1,class(fx));
                        if q.tgt(kn1,7)
                            dat.fyr{kn}(:,ks)=nanfill;
                        end
                        dat.vxr{kn}(:,ks)=nanfill;
                        dat.vyr{kn}(:,ks)=nanfill;
                        dat.pxr{kn}(:,ks)=nanfill;
                        dat.pyr{kn}(:,ks)=nanfill;
                        dat.ms{kn}(:,ks)=nanfill;
                        if includeft
                            dat.fy_ft{kn}(:,ks)=nanfill;
                            dat.fx_ft{kn}(:,ks)=nanfill;
                        end
                        %     dat.vxr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %     dat.vyr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %     dat.pxr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %     dat.pyr{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %     dat.ms{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %     if includeft
                        %         dat.fy_ft{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %         dat.fx_ft{kn}(:,ks)=nan(2*wsize+1,1,class(fx));
                        %     end
                        % end
                        
                        nanval = nan('like',fx);
                        dat.mtime(kn,ks) = nan;
                        dat.pathlen(kn,ks) = nanval;
                        dat.mdist(kn,ks) = nanval;
                        dat.pathlen_rel(kn,ks) = nanval;
                        dat.vtmax(kn,ks) = nanval;
                        dat.vxmax(kn,ks) = nanval;
                        dat.vymax(kn,ks) = nanval;
                        dat.ra(kn,ks) = nan;
                        dat.tvtmax(kn,ks) = nan;
                        dat.tvxmax(kn,ks) = nan;
                        dat.n(kn,ks) = nan;
                    end
                else disp(['skipped (empty) #', num2str(kn)])
                end
            else disp(['skipped (existence) #__', num2str([kn,kn1,sum(info.blocksize(1:kb-1))])]);  keyboard
            end
%             figure,
%             hold on
%             plot(15*vxr)
%             plot(fyr)
%             vline(tvxmax,'k','linewidth',2)
%             vline(tvxmax-wsize,'r','linewidth',2)
%             vline(tvxmax+wsize,'r','linewidth',2)
%             vline(tvxmax-120,'g','linewidth',2)
%             vline(tvxmax+120,'g','linewidth',2)
%             was_good = dat.good(kn,ks);
%             was_nan = sum(isnan(dat.fyr{kn}(:,ks)));
%             title(['kn:' num2str(kn) ', good: ' num2str(was_good) ', nans: ' num2str(was_nan)])
%             keyboard
            if 0
                figure,
                plot(ms)
            end
        end;
    end;
end

% f=fieldnames(dat);
return;


function [xr,yr]=rotate2(phi, x, y)
xr = cos(phi)*x - sin(phi)*y;
yr = sin(phi)*x + cos(phi)*y;
return
