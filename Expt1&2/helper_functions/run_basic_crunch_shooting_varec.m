function [dat,info]=run_basic_crunch_shooting_varec(dir1)
% I've altered this code to ignore the first movement, since it is a
% passive movement anyway to get subjects to the bottom starting position.
% It should only be used in such circumstances; otherwise you are losing
% the first data point from every block.

% for force and position etc. variability comparison

% window size of 225*2+1

dbstop if error
info.base_dir = dir1;
info.FFMAG=15;
tic;
d=dir([info.base_dir,'\*.mat'])

wsize=225;
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
info.blockletter = unique(filename_list_blockid);

dt=dir([info.base_dir,'\*.tgt']);
if length(dt) > 0,
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
    info.blocksize = mode(tgt_dat_length,1);
    tgt_all = vertcat(tgt_dat{1,:});
    for k2=1:length(info.blockletter),
        for k1a=1:length(info.sublist), for k1b=1:length(info.sublist),
                if length(tgt_dat{k1a,k2}(:)) == length(tgt_dat{k1b,k2}(:)),
                    tgt_dat_equal(k1a,k1b,k2) = all(tgt_dat{k1a,k2}(:) == tgt_dat{k1b,k2}(:));
                else tgt_dat_equal(k1a,k1b,k2) = 0;
                end
                
            end; end
    end
    tgt_dat_equal_summary = all(reshape(tgt_dat_equal,k1a*k1b,k2));
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

info

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
for kb = 1:length(info.blocksize), disp(['Block ',info.blockletter(kb)]); % block counter
    for ks = 1:length(info.sublist), disp([info.sublist{ks},' ',num2str(toc)]) % subject counter
        filename = [info.base_dir,'\',info.sublist{ks},info.blockletter{kb},'.mat'];
        q = load(filename);
        if 1, %~isfield(q,'tgt'),
            cs = [0,cumsum(info.blocksize)];
            q.tgt = tgt_all(cs(kb)+1:cs(kb+1),:);
        end
        if min(size(q.tgt))==1, q.tgt = reshape(q.tgt,8,[])'; end;
        if length(q.tgt)~=info.blocksize(kb),
            disp('unexpected tgt fle length');
        end
        if max(q.tgt(:,1))<=32, NUM_DIRS=32; else NUM_DIRS=360; end
        phi = -q.tgt(:,1)*2*pi/NUM_DIRS;
        x0_tgt = [0; cumsum(TARGET_DIST*cos(-phi))];
        y0_tgt = [0; cumsum(TARGET_DIST*sin(-phi))];
        % kn is overall movement number)
        for kn1=2:info.blocksize(kb), %The first 'movement' is just to get to the initial target.
            kn=kn1+sum(info.blocksize(1:kb-1)); %if ks == 2 && kn == 1072, keyboard, end
            if length(q.rawdata)>=kn1,
                if ~isempty(q.rawdata{kn1}) && size(q.rawdata{kn1},1)>(n1+10)
                    x = q.rawdata{kn1}(n1:end,1); y = q.rawdata{kn1}(n1:end,2);
                    vx = q.rawdata{kn1}(n1:end,3); vy = q.rawdata{kn1}(n1:end,4);
                    ax = q.rawdata{kn1}(n1:end,5); ay = q.rawdata{kn1}(n1:end,6);
                    fx = q.rawdata{kn1}(n1:end,9); fy = q.rawdata{kn1}(n1:end,10);
                    timing_data = q.rawdata{kn1}(n1:end,13);
                    
                    % buffer fx and fy
                    firstnonzero=find(fx~=0,1);
                    fx(1:firstnonzero)=nanmean(fx((firstnonzero+1):(firstnonzero+10)));
                    firstnonzero=find(fy~=0,1);
                    fy(1:firstnonzero)=nanmean(fy((firstnonzero+1):(firstnonzero+10)));
                    
                    fx_ft = q.rawdata{kn1}(n1:end,7); fy_ft = q.rawdata{kn1}(n1:end,8);
                    %[xr,yr,vxr,vyr] = rotate_data(q);
                    
                    if isfield(q,'tgt'),
                        [xr,yr] = rotate2(phi(kn1), x-x0_tgt(kn1), y-y0_tgt(kn1));
                        [vxr,vyr] = rotate2(phi(kn1), vx, vy);
                        [axr,ayr] = rotate2(phi(kn1), ax, ay);
                        [fxr,fyr] = rotate2(phi(kn1), fx, fy);
                        [fxr_ft,fyr_ft] = rotate2(phi(kn1), fx_ft, fy_ft);
                    else %assume 90/270
                        if rem(kn1,2)==1,
                            xr=-y; yr=x; vxr=-vy; vyr=vx; axr=-ay; ayr=ax; fxr=-fy; fyr=fx; %odd mvts (down)
                        else
                            xr=y+0.1; yr=-x; vxr=vy; vyr=-vx; axr=ay; ayr=-ax; fxr=fy; fyr=-fx; %even mvts (up)
                        end
                    end
                    if phi>(pi/2+.0001) & phi < (1.5*pi-0.001),
                        keyboard;
                    end
                    
                    % t40 is last time at which movement passes 4 cm
                    t40 = max(find([[xr(1:end-1)<=0.04]&[xr(2:end)>0.04]]));
                    if t40 >200,
                        vxr(1:t40-200)=nanmean(vxr((t40-200):(t40-190)));
                        vyr(1:t40-200)=nanmean(vyr((t40-200):(t40-190)));
                        % Andrew added the next two lines, Feb 2014.
                        xr(1:t40-200)=nan;
                        yr(1:t40-200)=nan;
                    end
                    
                    vt = sqrt(vxr.^2+vyr.^2);
                    % qra is number of samples where person is moving faster
                    % than threshold
                    qra = cumsum(vt>v_th); % v_th is threshold for movement (.05)
                    ra=min(find([qra(v_tht+1:end)-qra(1:end-v_tht)==v_tht]));
                    % ra is last time the person is under threshold
                    
                    % if person is begins too early or too late then it's not good
                    if isempty(ra) | ra>7000, dat.good(kn,ks) = 0; message='start late';
                    elseif ra<20, dat.good(kn,ks) = 0; message=['start early: ra=',num2str(ra)];
                    elseif isempty(t40), dat.good(kn,ks) = 0; message='miss 4cm';
                    else dat.good(kn,ks) = 1;
                    end
                    
                    % if movement is too fast or too slow (velthresh)
                    vmaxthresh = 1.6; %.5;
                    vminthresh=.23;
                    dat.timing_data{kn,ks}=timing_data;
                    if isempty(vxr)
                        dat.good(kn,ks) = 0;
                        message='movement too short';
                    elseif max(abs(vxr))>vmaxthresh | max(abs(vxr))<vminthresh
                        dat.vmax(kn,ks)=max(abs(vxr));
                        dat.good(kn,ks)=0;
                        if max(abs(vxr))>vmaxthresh, message=['too fast: ',num2str(max(abs(vxr)),3)];
                        else message=['too slow: ',num2str(max(abs(vxr)),2)];
                        end
                    else 
                        dat.vmax(kn,ks)=max(abs(vxr));
                    end
                    
                    %                    figure(99); subplot(211); plot([vxr,vyr]); title(num2str([kn,kn1,ks,ra,t40])); subplot(212); plot([xr,yr]); shg; pause;
                    if  dat.good(kn,ks),
                        
                        % t0 where the actual movement starts
                        t0 = max(1,ra-t_startoffset);
                        tvtmax = min(find(vt(t0:end)==max(vt(t0:end))))+t0-1; % time of first max vel
                        tvxmax = min(find(abs(vxr(t0:end))==max(abs(vxr(t0:end)))))+t0-1;
                        tvymax = min(find(abs(vyr(t0:end))==max(abs(vyr(t0:end)))))+t0-1;
                        
                        tf = min(length(vxr), t0+2*(tvtmax-t0));
                        % end of movement, either the total length or twice the
                        % time from t0 to tvtmax after t0
                        
                        dat.mtime(kn,ks) = tf-t0; %tstop-ra; %tf-t0;
                        dat.pathlen(kn,ks) = sum(vt)*5;
                        dat.mdist(kn,ks) = sqrt([xr(tf)-xr(t0)]^2+[yr(tf)-yr(t0)]^2);
                        dat.pathlen_rel(kn,ks) = dat.pathlen(kn,ks)/dat.mdist(kn,ks);
                        dat.vtmax(kn,ks) = max(vt);
                        dat.vxmax(kn,ks) = vxr(tvxmax);
                        dat.vymax(kn,ks) = vyr(tvymax);
                        dat.ra(kn,ks) = ra;
                        dat.tvtmax(kn,ks) = tvtmax - ra;
                        dat.tvxmax(kn,ks) = tvxmax - ra;
                        dat.n(kn,ks) = length(vxr);
                        if q.tgt(kn1,7),               %if FC trial then save force data
                            if tvxmax > wsize, % time of max x velocity is more than windowsize
                                if tvxmax+wsize > length(fyr),
                                    Nmissing = tvxmax+wsize - length(fyr);
                                    
                                    dat.fyr{kn}(:,ks)=nan(2*wsize+1,1);
                                    dat.vxr{kn}(:,ks)=nan(2*wsize+1,1);
                                    dat.vyr{kn}(:,ks)=nan(2*wsize+1,1);
                                    dat.pxr{kn}(:,ks)=nan(2*wsize+1,1);
                                    dat.pyr{kn}(:,ks)=nan(2*wsize+1,1);
                                    
                                    dat.fyr{kn}(1:wsize*2+1-Nmissing,ks)=fyr(tvxmax-wsize:end);
                                    dat.vxr{kn}(1:wsize*2+1-Nmissing,ks)=vxr(tvxmax-wsize:end);
                                    dat.vyr{kn}(1:wsize*2+1-Nmissing,ks)=vyr(tvxmax-wsize:end);
                                    dat.pxr{kn}(1:wsize*2+1-Nmissing,ks)=xr(tvxmax-wsize:end);
                                    dat.pyr{kn}(1:wsize*2+1-Nmissing,ks)=yr(tvxmax-wsize:end);
                                else
                                    dat.fyr{kn}(:,ks) = fyr(tvxmax + [-wsize:wsize]);
                                    dat.vxr{kn}(:,ks) = vxr(tvxmax + [-wsize:wsize]);
                                    dat.vyr{kn}(:,ks) = vyr(tvxmax + [-wsize:wsize]);
                                    
                                    %                            if mod(kn,2)==1
                                    dat.pxr{kn}(:,ks) = xr(tvxmax + [-wsize:wsize]);
                                    dat.pyr{kn}(:,ks) = yr(tvxmax + [-wsize:wsize]);
                                    %                            else
                                    %                                dat.pxr{kn}(:,ks) = -xr(tvxmax + [-wsize:wsize])+.1;
                                    %                            end
                                end
                            else
                                if tvxmax+wsize > length(fyr)
                                    Ndata = length(fyr);
                                    dat.fyr{kn}(:,ks)=nan(2*wsize+1,1);
                                    dat.vxr{kn}(:,ks)=nan(2*wsize+1,1);
                                    dat.pxr{kn}(:,ks)=nan(2*wsize+1,1);
                                    dat.vyr{kn}(:,ks)=nan(2*wsize+1,1);
                                    dat.pyr{kn}(:,ks)=nan(2*wsize+1,1);
                                    
                                    dat.fyr{kn}(wsize+2-tvxmax:wsize+1-tvxmax+Ndata,ks)=fyr;
                                    dat.vxr{kn}(wsize+2-tvxmax:wsize+1-tvxmax+Ndata,ks)=vxr;
                                    dat.pxr{kn}(wsize+2-tvxmax:wsize+1-tvxmax+Ndata,ks)=xr;
                                    dat.vyr{kn}(wsize+2-tvxmax:wsize+1-tvxmax+Ndata,ks)=vyr;
                                    dat.pyr{kn}(wsize+2-tvxmax:wsize+1-tvxmax+Ndata,ks)=yr;
                                    
                                else
                                    % pad vxr and fyr in front
                                    padsize=wsize-tvxmax+1;
                                    padvalfyr=nan;%nanmean(fyr(find(fyr~=0,10)));
                                    padvalvxr=nan;%nanmean(vxr(find(vxr~=0,10)));
                                    padvalvyr=nan;%nanmean(vyr(find(vyr~=0,10)));
                                    %                                if mod(kn,2)==1
                                    padvalxr=nan;%nanmean(xr(find(xr~=0,10)));
                                    padvalyr=nan;%nanmean(yr(find(yr~=0,10)));
                                    %                                else
                                    %                                    padvalxr=-nanmean(xr(1:10))+.1;
                                    %                                end
                                    dat.fyr{kn}(1:padsize,ks)=padvalfyr*ones(padsize,1);
                                    dat.vxr{kn}(1:padsize,ks)=padvalvxr*ones(padsize,1);
                                    dat.pxr{kn}(1:padsize,ks)=padvalxr*ones(padsize,1);
                                    dat.vyr{kn}(1:padsize,ks)=padvalvyr*ones(padsize,1);
                                    dat.pyr{kn}(1:padsize,ks)=padvalyr*ones(padsize,1);
                                    
                                    dat.fyr{kn}((padsize+1):(wsize*2+1),ks) = fyr(1:(tvxmax + wsize));
                                    dat.vxr{kn}((padsize+1):(wsize*2+1),ks) = vxr(1:(tvxmax + wsize));
                                    dat.vyr{kn}((padsize+1):(wsize*2+1),ks) = vyr(1:(tvxmax + wsize));
                                    
                                    %                                if mod(kn,2)==1
                                    dat.pxr{kn}(padsize+1:wsize*2+1,ks) = xr(1:tvxmax + wsize);
                                    dat.pyr{kn}(padsize+1:wsize*2+1,ks) = yr(1:tvxmax + wsize);
                                    %                                else
                                    %                                    dat.pxr{kn}(padsize+1:wsize*2+1,ks) = -xr(1:tvxmax + wsize)+.1;
                                    %                                end
                                end
                            end
                        else
                            % save velocity and position data if non error clamp
                            if tvxmax > wsize, % time of max x velocity is more than windowsize
                                if tvxmax+wsize > length(vxr)
                                    Nmissing = tvxmax+wsize - length(vxr);
                                    dat.vxr{kn}(:,ks)=nan(2*wsize+1,1);
                                    dat.vyr{kn}(:,ks)=nan(2*wsize+1,1);
                                    dat.pxr{kn}(:,ks)=nan(2*wsize+1,1);
                                    dat.pyr{kn}(:,ks)=nan(2*wsize+1,1);

                                    dat.vxr{kn}(1:wsize*2+1-Nmissing,ks)=vxr(tvxmax-wsize:end);
                                    dat.vyr{kn}(1:wsize*2+1-Nmissing,ks)=vyr(tvxmax-wsize:end);
                                    dat.pxr{kn}(1:wsize*2+1-Nmissing,ks)=xr(tvxmax-wsize:end);
                                    dat.pyr{kn}(1:wsize*2+1-Nmissing,ks)=yr(tvxmax-wsize:end);
                                else
                                    dat.vxr{kn}(:,ks) = vxr(tvxmax + [-wsize:wsize]);
                                    dat.vyr{kn}(:,ks) = vyr(tvxmax + [-wsize:wsize]);
                                    
                                    %                            if mod(kn,2)==1
                                    dat.pxr{kn}(:,ks) = xr(tvxmax + [-wsize:wsize]);
                                    dat.pyr{kn}(:,ks) = yr(tvxmax + [-wsize:wsize]);
                                    %                            else
                                    %                                dat.pxr{kn}(:,ks) = -xr(tvxmax + [-wsize:wsize])+.1;
                                    %                            end
                                end
                            else
                                if tvxmax+wsize > length(vxr)
                                    Ndata = length(vxr);
                                    dat.vxr{kn}(:,ks)=nan(2*wsize+1,1);
                                    dat.pxr{kn}(:,ks)=nan(2*wsize+1,1);
                                    dat.vyr{kn}(:,ks)=nan(2*wsize+1,1);
                                    dat.pyr{kn}(:,ks)=nan(2*wsize+1,1);
                                    
                                    dat.vxr{kn}(wsize+2-tvxmax:wsize+1-tvxmax+Ndata,ks)=vxr;
                                    dat.pxr{kn}(wsize+2-tvxmax:wsize+1-tvxmax+Ndata,ks)=xr;
                                    dat.vyr{kn}(wsize+2-tvxmax:wsize+1-tvxmax+Ndata,ks)=vyr;
                                    dat.pyr{kn}(wsize+2-tvxmax:wsize+1-tvxmax+Ndata,ks)=yr;
                                else
                                    % pad vxr and fyr in front
                                    padsize=wsize-tvxmax+1;
                                    padvalvxr=nan;%nanmean(vxr(find(vxr~=0,10)));
                                    padvalvyr=nan;%nanmean(vyr(find(vyr~=0,10)));
                                    %                                if mod(kn,2)==1
                                    padvalxr=nan;%nanmean(xr(find(xr~=0,10)));
                                    padvalyr=nan;%nanmean(yr(find(yr~=0,10)));
                                    %                                else
                                    %                                    padvalxr=-nanmean(xr(1:10))+.1;
                                    %                                end
                                    dat.vxr{kn}(1:padsize,ks)=padvalvxr*ones(padsize,1);
                                    dat.pxr{kn}(1:padsize,ks)=padvalxr*ones(padsize,1);
                                    dat.vyr{kn}(1:padsize,ks)=padvalvyr*ones(padsize,1);
                                    dat.pyr{kn}(1:padsize,ks)=padvalyr*ones(padsize,1);
                                    
                                    dat.vxr{kn}((padsize+1):(wsize*2+1),ks) = vxr(1:(tvxmax + wsize));
                                    dat.vyr{kn}((padsize+1):(wsize*2+1),ks) = vyr(1:(tvxmax + wsize));
                                    
                                    %                                if mod(kn,2)==1
                                    dat.pxr{kn}(padsize+1:wsize*2+1,ks) = xr(1:tvxmax + wsize);
                                    dat.pyr{kn}(padsize+1:wsize*2+1,ks) = yr(1:tvxmax + wsize);
                                    %                                else
                                    %                                    dat.pxr{kn}(padsize+1:wsize*2+1,ks) = -xr(1:tvxmax + wsize)+.1;
                                    %                                end
                                end
                            end
                        end
                    else
                        if mod(kn,2)==0, disp(['aborted (badtrial) #', num2str(kn),': ',message]), end
                        if q.tgt(kn1,7)
                        dat.fyr{kn}(:,ks)=nan(2*wsize+1,1);
                        dat.vxr{kn}(:,ks)=nan(2*wsize+1,1);
                        dat.vyr{kn}(:,ks)=nan(2*wsize+1,1);
                        dat.pxr{kn}(:,ks)=nan(2*wsize+1,1);
                        dat.pyr{kn}(:,ks)=nan(2*wsize+1,1);
                        else
                        dat.vxr{kn}(:,ks)=nan(2*wsize+1,1);
                        dat.vyr{kn}(:,ks)=nan(2*wsize+1,1);
                        dat.pxr{kn}(:,ks)=nan(2*wsize+1,1);
                        dat.pyr{kn}(:,ks)=nan(2*wsize+1,1);    
                        end
                    end
                else
                    if mod(kn,2)==0, disp(['skipped (empty) #', num2str(kn)]), end
                end
            else disp(['skipped (existence) #__', num2str([kn,kn1,sum(info.blocksize(1:kb-1))])]);  keyboard
            end
        end;
    end;
end

f=fieldnames(dat);
return;


function [xr,yr]=rotate2(phi, x, y);
xr = cos(phi)*x - sin(phi)*y;
yr = sin(phi)*x + cos(phi)*y;
return
