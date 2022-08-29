function matfilename = convert_imt_gz2mat2(gzfilename,gzfiledir,matfiledir, matfile_dataprecision) 
%
% matfilename = convert_imt_gz2mat2(gzfilename,gzfiledir,matfiledir);
%
% converts gz files containing one particpant's data into matlab mat files
% named something like: 'xxxa.mat' (one mat file per set)  
% it assumes that the first 3 chars in each gz filename are the subject id
% data files for each are put into a cell array named 'rawdata'
% if available, the tgt file is put into an array named 'tgt'
% conversion should take about a minute per participant
% ........
% inputs/arguments:
% gzfilename: name of gzip file to convert, e.g. z_xxx.gz
%             '*' --> all gz files converted
% gzfiledir: directory where gz file resides
% matfiledir: directory for mat files
%     note: these 1st 3 inputs are all string variables 
% matfile_dataprecision: single (32 bit) or double precision (64 bit) for mat-file data
%                        1 --> single (default), 2 --> double 
% ........
% defaults
% if not specified gzfilename defaults to '*' (all gz files in directory)
% if not specified gzfiledir defaults to the current directory
% if not specified matfiledir defaults to gzfiledir
% if not specified matfile_dataprecision defaults to 1
% ........
% MAS 4/17/09
%
% version 2:  
% (1) uses matlab's new internal support for gz and tar archives versus an external program like alzip
% (2) provides an option for the precision of matlab mat-files - and uses single presicion (32-bit) 
%      by default to save disk space because double precision (64-bit) data types (matlab default) 
%      are overkill for our data given that the robot's sensors / ADC are intrinsically only 16 bits. 
%
% % to do:  
% (1) remove pesky for-loop index warning:
% Warning: FOR loop index with character indexing will return character
%          instead of double in future release.
% (2) send progress display to a figure window?? can we have a scrolling dsiplay??
% (3) allow for updating directories --> i.e. only process new gz files (ones without mat files)

if nargin < 1, gzfilename='*'; end
if nargin < 2, gzfiledir=pwd; end
if nargin < 3, matfiledir=gzfiledir; end
if isempty(matfiledir), matfiledir=gzfiledir; end
if nargin < 4, matfile_dataprecision = 1; end

cd(gzfiledir)
if gzfilename(1)=='*', gzfilenamelist = cellstr(findfile('*','dir *.gz')); 
else gzfilenamelist={gzfilename};
end

tic;
for k=1:length(gzfilenamelist),
    gzfilename = gzfilenamelist{k};
    gzfileprefix = gzfilename(1:max(find(gzfilename=='.'))-1);
%         keyboard
    subject_id = gzfilename(1:3);
    tempdir = [gzfiledir,'\temp'];
    tempdir1 = [tempdir,'\',gzfileprefix];

    if exist(tempdir1)==7, 
        cd(tempdir1);
        delete('*.dat') %clean-up temp dir
        delete('*.tar')
        delete('*.asv')
        delete('*')
        cd ..
        rmdir(tempdir1);
    end
    if exist(tempdir)==7
        cd(tempdir)
        delete('*.dat') %clean-up temp dir
        delete('*.tar')
        delete('*.asv')
        delete('*')
        cd ..
        rmdir(tempdir);
    end
    cd(gzfiledir)
    gunzip(gzfilename,tempdir);
    cd(tempdir)
%     keyboard
    
    untar(gzfileprefix,'.\temp');
    cd('temp')
    newfile = findfile('*','dir *a1.dat');
    newfile1 = deblank(newfile(1,:));
    prefix  = newfile1(1,1:end-6);                 %look at this later
    for ks=1:26,  
        tt.a=toc;
        s = char('a'+ks-1);  %  1:26 --> a:z 
        maxfilenum = length(dir([prefix,char(s),'*.dat']));
        if maxfilenum>0,
            clear rawdata tgt;
            tt.b=toc;
            for k1=1:maxfilenum, 
                datfile=[prefix,char(s),num2str(k1),'.dat'];
                if exist(datfile), 
                    if  matfile_dataprecision == 1, rawdata{k1}=single(load(datfile)); end
                    if  matfile_dataprecision == 2, rawdata{k1}=load(datfile); end
                else
                    disp(['missing file: ',datfile]);
                end; 
            end
            tt.c=toc;
            tgtfile = findfile('*',['dir *',s,'.tgt']);
            matfile = [matfiledir,'\',subject_id,s,'.mat'];
            if (size(tgtfile,1)==1), 
                tgt = load(tgtfile);
                save(matfile,'rawdata','tgt'); 
            elseif (size(tgtfile,1)>1),
                tgt = load(tgtfile(1,:));
                save(matfile,'rawdata','tgt'); 
            elseif (size(tgtfile,1)==0),
                save(matfile,'rawdata'); 
            end
            %tt.d = toc; disp([[tt.a, tt.b, tt.c, tt.d], diff(fliplr([tt.a, tt.b, tt.c, tt.d]))]);
            disp(['matfile name: ',matfile, '   # of dat files: ',num2str(maxfilenum), '   # of tgt files: ',num2str(size(tgtfile,1))]);
        end
    end
    delete('*.dat') %clean-up temp dir
    delete('*.tgt')
    delete('*.tar')
    delete('*.asv')
    delete('*')
    length(findfile('*','dir *.dat'));
    cd ..
    rmdir('temp');
    delete(gzfileprefix)
    cd ..
    try
        rmdir('temp')
    catch
        display('rmdir second try')
        try
            rmdir('temp')
        catch
            display('rmdir second fail')
            keyboard
        end
    end
    matfilename =[subject_id,s,'.mat'];
    cd ..
end
return

function [names, sizes] = findfile(comp_str, dir_cmd)

if nargin < 2, dir_cmd = 'dir'; end
[s,w] = dos([dir_cmd, '/-c']);
nw = find(w==10); % find newline characters in return string
nwb = nw(5:(length(nw)-3));
nwe = nw(6:(length(nw)-2));
clear names sizes;
if length(nwb)>0,
    for k = 1:length(nwb),   
        names(k,1:(nwe(k)-nwb(k)-40)) = w((nwb(k)+40):(nwe(k)-1));
        sizes(k) = max([ 0, str2num(w([(nwb(k)+25):(nwb(k)+38)])) ]);
    end
    if ( deblank(names(1,:)) == '.' ),
        names = names(3:size(names,1),:);
    end   
    if nargin > 0;
        name_ind = ones(size(names,1),1);
        for k = 1:length(comp_str),
            if ismember(comp_str(k),'#*'),
                if comp_str(k) == '#',
                    name_ind = name_ind .* ismember(names(:,k),'0123456789');
                end   
            else name_ind = name_ind .* ismember(lower(names(:,k)),lower(comp_str(k)));
            end
        end
        names = lower(names(find(name_ind),:));   
        sizes = sizes(find(name_ind));   
    end
else
    names=[]; sizes=[];
end
