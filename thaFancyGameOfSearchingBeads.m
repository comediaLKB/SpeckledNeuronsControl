%% %%%%%%%%%%%%%%%%%%
% assuming having the map "TT" from "calibrate_beads_sample.m"
% the fancy game of sreaching beads
if ~exist('TT')
    load('D:\Users\Comedia\moro\06052019\TT.mat')
end
%%
% check if the cameras has been loaded yet
% delete(PCOcam);
% delete(BASLERcam);
% if ~exist('PCOcam'); end
% if ~exist('BASLERcam'); ; end
% PCOcam = PCOcamera(200);
% BASLERcam = BASLERcamera(200);
%% epifluiorescence from below
BASLERexp = 0.1;
fov_size = 1;
intensity = 0.5;

BASLERsrc.ExposureTime = BASLERexp*1000;

fig = figure('position', [200, 200, 600, 600]);

pause(0.1);
illumination = drawpoints(DMD_SIZE, DMD_SIZE/2, min(DMD_SIZE)*fov_size)*intensity;
blankscreen = rot90(zeros(DMD_SIZE));

while ishandle(fig),
    fullscreen(illumination, 2);
    orbeads = getsnapshot(BASLERvid);
%     pause(0.3);
    fullscreen(blankscreen,2);
    imagesc(orbeads);
    pbaspect([1 1 1]);
    drawnow
end



%% SEARCH BEADS
PCOexp = 2000;
BASLERexp = 0.05;

BASLERsrc.ExposureTime = BASLERexp*1000;
PCOsrc.E2ExposureTime = PCOexp;

disks_d = 20;
% FOV in percentage of the DMD
fov_size = 0.8;

% [XX,YY] = ndgrid(1:DMD_SIZE(1),1:DMD_SIZE(2));
% illumination = (XX-DMD_SIZE(1)/2).^2 + (YY-DMD_SIZE(2)/2).^2 <= (DMD_SIZE(2)*fov*0.5)^2;
% 

fig = figure('position', [200, 200, 1200, 600]);
pause(0.1);
se = strel('disk',5);
illumination = drawpoints(DMD_SIZE, DMD_SIZE/2, min(DMD_SIZE)*fov_size);
blankscreen = rot90(zeros(DMD_SIZE));

while ishandle(fig),
    % acquire
    fullscreen(illumination, 2);
    orbeads = getsnapshot(BASLERvid);  
    topbeads = getsnapshot(PCOvid);
    fullscreen(blankscreen,2);

    subplot(1,2,1);
    imagesc(orbeads);
    pbaspect([1 1 1]);
    subplot(1,2,2);
    imagesc(topbeads);
    pbaspect([1 1 1]);
%     refresh;
    drawnow
end



%% SEARCH BEADS FINDING THE CENTERS
% PCOexp = 100;
% BASLERexp = 0.1;
% 
% BASLERsrc.ExposureTime = BASLERexp*1000;
% PCOsrc.E2ExposureTime = PCOexp;
% 
dmin = 5;
dmax = 15;

% % FOV in percentage of the DMD
% fov_size = 0.8;

% [XX,YY] = ndgrid(1:DMD_SIZE(1),1:DMD_SIZE(2));
% illumination = (XX-DMD_SIZE(1)/2).^2 + (YY-DMD_SIZE(2)/2).^2 <= (DMD_SIZE(2)*fov*0.5)^2;
% 

fig = figure('position', [200, 200, 600, 600]);
% % bottom camera exposure
% srcbot.ExposureTime = 1000;
% % sCMOS exposure
% exposure = 200; % max 2s
% multiple_frames_sum = 1;
% camera.setMultipleFrames(multiple_frames_sum);
% camera.setExposure(exposure);



pause(0.1);
% se = strel('disk',disks_d);
illumination = drawpoints(DMD_SIZE, DMD_SIZE/2, min(DMD_SIZE)*fov_size);
blankscreen = rot90(zeros(DMD_SIZE));

while ishandle(fig),
    % acquire
%     fullscreen(rot90(illumination),2);
    fullscreen(illumination, 2);
    orbeads = getsnapshot(BASLERvid);
%     orbeads = BASLERcam.acquireOne;    
    topbeads = getsnapshot(PCOvid);
%     topbeads = PCOcam.acquireOne;
    fullscreen(blankscreen,2);
    % 40 px diameter the beads, more or less
    % try to find the circles woth imfindcircles
    % visualize
    subplot(2,2,1);
    imagesc(orbeads);
    pbaspect([1 1 1]);
    subplot(2,2,2);
    imagesc(orbeads);
    
    [centers, radii, metric] = imfindcircles(orbeads,[dmin dmax]);
    viscircles(centers, radii,'EdgeColor','r');
    text1 = text(20,40,'','Color','white');
    set(text1,'String',sprintf('%d beads',size(centers, 1)));
    subplot(2,2,3);
    imagesc(topbeads);
%     refresh;
    drawnow
end
profile off
beads_pos = [];
for pt = centers'
    v = TT*[pt(1) pt(2) 1]';
    beads_pos = [beads_pos [v(1),v(2)]'];
end
beads_pos = beads_pos';
% record the original beads position for future analisys (to be reloaded in
% case)
original_beads_pos = beads_pos;
%%
beads_pos = original_beads_pos;
%% show the beads position with number
% show the position in a screen with the numbers, so then one can 
% remove selectively
figure();
tmpi = drawpoints(DMD_SIZE, beads_pos, disks_d);
imshow(tmpi(1:3:end,1:3:end));
for idx = 1:size(beads_pos,1)
    pt = beads_pos(idx,:);
    text(pt(1)/3+4, pt(2)/3, int2str(idx), 'color', 'green')
end
%% delete some beads
todelete = [21 22];
beads_pos(todelete,:) = [];
disp(length(beads_pos));

%% select a subset of beads as background
fg_beads_no = 20;
% bg_beads_no = length(beads_pos)-fg_beads_no;
bg_beads_no = 40;
iids = randperm(length(beads_pos));
iids = iids(1:bg_beads_no);

bg_beads_pos = beads_pos(iids,:);
beads_pos(iids,:) = [];
beads_pos(fg_beads_no+1:end,:) = [];


% show the beads position and the bg
figure();
tmpi = drawpoints(DMD_SIZE, beads_pos, disks_d) + drawpoints(DMD_SIZE, bg_beads_pos, disks_d)*0.5;
imshow(tmpi(1:3:end,1:3:end));

%% select a subset of random beads checking if they emit something detectable
% select a ROI
disk_diam_px = 20;
PCOexp = 100;
PCOsrc.E2ExposureTime = PCOexp;
fullscreen(drawpoints(DMD_SIZE, beads_pos, disk_diam_px), 2);
data = getsnapshot(PCOvid);
pause(PCOexp/1000);
fullscreen(rot90(zeros(DMD_SIZE)),2);
figure(10);imagesc(data);
h = imrect;
crop = round(h.getPosition);
frame_size = [crop(3),crop(4)];
crop(3:4) = crop(3:4)+crop(1:2);
imagesc(data(crop(2):crop(4),crop(1):crop(3)));
%% measure a sample of some beads and detect the threshold observing the
% distribution, then take some random ones wich are above the noise
all_bs = length(beads_pos);
iids = randperm(length(beads_pos));
contrast = [];
variance = [];
for ids = iids(1:all_bs)
    fullscreen(drawpoints(DMD_SIZE, beads_pos(ids,:), disk_diam_px), 2);
    data = getsnapshot(PCOvid);
    pause(PCOexp/1000);
    data = data(crop(2):crop(4),crop(1):crop(3));
    fullscreen(rot90(zeros(DMD_SIZE)),2);
    % contrast
    contrast = [contrast std(single(data(:)))/mean(data(:))];
    % variance
    variance = [variance std(single(data(:)))^2];
end

% find the subset of beads in the central part of the distribution

% pickup 20 of these index, and create the new beads_pos and bg_beads_pos



%% acquire a sequence in tiff file
% with a tiff file fot the top image
% a tiff file with the bottom image (for visualization purposes)
% a dump of the points in the DMD space

PCOexp = 1500;
BASLERsrc.ExposureTime = BASLERexp*1000;
PCOsrc.E2ExposureTime = PCOexp;

randseed = int2str(rand()*10000);
top_tiff_file = strcat(root_folder,'video_',datestr(now,'ddmmyyyy'),'_', randseed, '_speckle.tiff');
bot_tiff_file = strcat(root_folder,'video_',datestr(now,'ddmmyyyy'),'_', randseed, '_bottomcamera.tiff');
point_csv_file = strcat(root_folder,'video_',datestr(now,'ddmmyyyy'),'_', randseed, '_spots.csv');

% exposure = 1000; % max 2s
% multiple_frames_sum = 2;
% camera.setMultipleFrames(multiple_frames_sum);
% camera.setExposure(exposure);
% 
% srcbot.ExposureTime = 10000;

previewfig= figure('Name','preview','NumberTitle','off', 'position', [300, 300, 800, 400]);

for pt = beads_pos'
    fullscreen(drawpoints(DMD_SIZE, pt, 20), 2);
    
    pause(0.1);
    botdata = getsnapshot(BASLERvid);
    imwrite(botdata, bot_tiff_file, 'writemode', 'append')
    topdata = rot90(getsnapshot(PCOvid), 2);
    imwrite(topdata, top_tiff_file, 'writemode', 'append')  
    
    subplot(1,2,1);
    imagesc(botdata);
    pbaspect([1 1 1]);
    subplot(1,2,2);
    imagesc(topdata);
    pbaspect([1 1 1]);
    drawnow
end
fullscreen(rot90(zeros(DMD_SIZE)),2);
close(previewfig);
csvwrite(point_csv_file, beads_pos);

%% aquire a sequence with a random choice of 'ss' permutation of different sources illuminated together
PCOexp = 1600;
multiple_shots = 20;
% realizations for every number of beads illuminated simultaneously
ss = 10;

BASLERsrc.ExposureTime = BASLERexp*1000;
PCOsrc.E2ExposureTime = PCOexp;

randseed = int2str(rand()*10000);
% top_tiff_file = strcat(root_folder,'video_',datestr(now,'ddmmyyyy'),'_', randseed, '_speckle.tiff');
% bot_tiff_file = strcat(root_folder,'video_',datestr(now,'ddmmyyyy'),'_', randseed, '_bottomcamera.tiff');
top_mat_file = strcat(root_folder,'video_',datestr(now,'ddmmyyyy'),'_', randseed, '_speckle.mat');
bot_mat_file = strcat(root_folder,'video_',datestr(now,'ddmmyyyy'),'_', randseed, '_bottomcamera.mat');
pointno_mat_file = strcat(root_folder,'video_',datestr(now,'ddmmyyyy'),'_', randseed, '_spotsno.mat');
pointpos_mat_file = strcat(root_folder,'video_',datestr(now,'ddmmyyyy'),'_', randseed, '_spotspos.mat');


% exposure = 1000; % max 2s
% multiple_frames_sum = 2;
% camera.setMultipleFrames(multiple_frames_sum);
% camera.setExposure(exposure);
% 
% srcbot.ExposureTime = 10000;

previewfig= figure('Name','preview','NumberTitle','off', 'position', [300, 300, 800, 400]);

% go from 1 to the max number of beads, and every step choose 'n' beads to
% illuminate simultaneously. take 'ss' number of realization for every
% number

video_data_top = [];
video_data_bot = [];
ptno_data = {};

ptno_data_idx = 1;
for beadssim = 1:length(beads_pos)
    for idx = 1:ss
        indices = randperm(length(beads_pos));
        indices = indices(1:beadssim);
        ptts = beads_pos(indices,:);
        
        fullscreen(drawpoints(DMD_SIZE, ptts, 20), 2);

        pause(0.1);
        botdata = getsnapshot(BASLERvid);
        botdata = botdata(1:5:end, 1:5:end);
%         imwrite(botdata, bot_tiff_file, 'writemode', 'append')
        video_data_bot = cat(3, video_data_bot, botdata);
        topdata = rot90(getsnapshot(PCOvid), 2);
        if multiple_shots>1
            for iiii = 2:multiple_shots
                topdata = (topdata + rot90(getsnapshot(PCOvid), 2))/2;
            end
        end
        
%         imwrite(topdata, top_tiff_file, 'writemode', 'append') 
        video_data_top = cat(3, video_data_top, topdata);
%         ptno_data = [ptno_data beadssim];
%         ptno_data = cat(3, ptno_data, ptts);
        ptno_data{ptno_data_idx} = ptts;
        ptno_data_idx = ptno_data_idx+1;
        subplot(1,2,1);
        imagesc(botdata);
        pbaspect([1 1 1]);
        subplot(1,2,2);
        imagesc(topdata);
        pbaspect([1 1 1]);
        drawnow
    end
end
fullscreen(rot90(zeros(DMD_SIZE)),2);
close(previewfig);


save(top_mat_file,'video_data_top','-v7.3');
save(pointpos_mat_file, 'beads_pos','-v7.3');
save(pointno_mat_file, 'ptno_data','-v7.3');
save(bot_mat_file,'video_data_bot','-v7.3');

%% live update to check the spread of the speckle patterns using one of the selected beads
PCOexp = 2000;
PCOsrc.E2ExposureTime = PCOexp;

beadno = 5;
fullscreen(drawpoints(DMD_SIZE, beads_pos(beadno,:), 40), 2);

previewfig= figure('Name','preview','NumberTitle','off', 'position', [300, 300, 800, 400]);
index = 1;
timescale = 0;
meanfluo = 0;
tic;
while ishandle(previewfig),
    data = getsnapshot(PCOvid);
    subplot(121);
    imagesc(rot90(data, 2));
    pbaspect([1 1 1]);
    subplot(122);
    timescale(index) = toc;
    meanfluo(index) = mean(data(:));
    plot(timescale, meanfluo);
    drawnow
    index = index +1;
end

fullscreen(rot90(zeros(DMD_SIZE)),2);

%% see what are the different speckle patterns

previewfig= figure('Name','preview','NumberTitle','off', 'position', [300, 300, 800, 400]);

for pt = beads_pos'
    fullscreen(drawpoints(DMD_SIZE, pt, 20), 2);
    pause(0.1);  
    data = getsnapshot(BASLERvid);
    topdata = rot90(getsnapshot(PCOvid), 2);
    % plot
    subplot(1,2,1);
    imagesc(data);
    pbaspect([1 1 1]);
    subplot(1,2,2);
    imagesc(topdata);
    pbaspect([1 1 1]);
    drawnow
end

fullscreen(rot90(zeros(DMD_SIZE)),2);
