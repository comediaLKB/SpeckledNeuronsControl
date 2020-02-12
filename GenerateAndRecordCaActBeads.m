%%
% make all DMD functions available to this code
% includere il servo shutter. cosi' da evitare bleached

% make a "nice" interface?
% where all the procedures are in a separate button, parameters in separate
% boxes, and the acquisition and others data displayed while acquiring...


% tmp_folder = 'D:\tmp\';
root_folder = strcat('D:\Users\Comedia\moro\',datestr(now,'ddmmyyyy'),'\');
[status, msg, msgID] = mkdir(root_folder);
disp(msg);

addpath('D:\Users\Comedia\moro\git_code\SpeckledFluoImagingControl\');
addpath('D:\Users\Comedia\moro\git_code\SpeckledFluoImagingControl\DMDConnect-master');


DMD_SIZE = [1920,1080];


%% Initialize the DMD
dmd1 = DMD();
dmd1.setMode(0);
fullscreen(rot90(zeros(DMD_SIZE)),1);
pause(0.1);
closescreen();
pause(0.1);
fullscreen(rot90(zeros(DMD_SIZE)),2);
% dmd1.display(rot90(zeros(DMD_SIZE)));

%% initialize the two servos
% pause(2);
% stage_x = servo_thorlabs(83843398);
% stage_y = servo_thorlabs(83843405);

%% initialize cameras
imaqreset;

PCOexp = 400;
PCOmultipleframes = 1;
BASLERexp = 200;

% BASLERvid = videoinput('gentl', 1, 'Mono8'); % basler bw (ac2040-55um)
BASLERvid = videoinput('gentl', 1, 'BayerBG8');  % basler a colori
BASLERsrc = getselectedsource(BASLERvid);
BASLERvid.FramesPerTrigger = 1;
BASLERsrc.ExposureTime = BASLERexp*1000;

PCOvid = videoinput('pcocameraadaptor', 0, 'CameraLink');
PCOsrc = getselectedsource(PCOvid);
PCOvid.FramesPerTrigger = 1;
PCOsrc.E1ExposureTime_unit = 'ms';
PCOsrc.E2ExposureTime = PCOexp;

% PCOcam = PCOcamera(PCOexp);
% BASLERcam = BASLERcamera(BASLERexp);
%% align upper and lower path

disks_d = 15;
% FOV in percentage of the DMD
fov_size = 0.8;
se = strel('disk',5);
illumination = drawpoints(DMD_SIZE, DMD_SIZE/2, min(DMD_SIZE)*fov_size)*0.03;
fullscreen(illumination, 2);
%% black screem
blankscreen = rot90(zeros(DMD_SIZE));
fullscreen(blankscreen,2);

%% preview top
PCOexp = 400;
% PCOexp = 2000;
PCOsrc.E2ExposureTime = PCOexp;

PCOvid.FramesPerTrigger = 1;
preview(PCOvid);
%% take a snapshot
topimage = getsnapshot(PCOvid);
figure(); imshow(topimage);

%% preview bottom
BASLERvid.FramesPerTrigger = 1;
preview(BASLERvid);


%% GRID PATTERN
% check focus with grid pattern (if needed)
sq = 120;
dim = 0.05;
% fullscreen(rot90(checkerboard(sq/2, DMD_SIZE(1)/sq, DMD_SIZE(2)/sq)>0.5)*dim,2);
fullscreen(rot90(zeros(DMD_SIZE)),2);

previewfig= figure('Name','preview','NumberTitle','off', 'position', [300, 300, 800, 400]);

while ishandle(previewfig),
    imagesc(rot90(getsnapshot(PCOvid), 2));
    pbaspect([1 1 1]);
    drawnow
end

fullscreen(rot90(zeros(DMD_SIZE)),2);

%%
%
% HERE THE CALIBRATION and THE FUNNY GAME BUT BE DONE

%% CHECK SPECKLES
% check the speckles with the first of the points in projectedbeads
topcamera.setExposure(exposure);
fullscreen(projectedbeads==1, 2);
% topcamera.preview
figure(666);
% data = topcamera.acquireOne;
data = getsnapshot(PCOvid);
imagesc(data);
disp(max(data(:)));

%% Black and White
% 
fullscreen(rot90(ones(DMD_SIZE)),2);
fullscreen(rot90(zeros(DMD_SIZE)),2);
%% generate the activity

% parameters
frames = 500;

radius_px = 100;
min_dist_px = 10;
% instead here set a radius from the center of the FOV, not a square pad
disk_diam_px = 20;
time_stepsize = 0.1;
% should be placed an average spikes no, or better an average 
% per second, and then a randn or similar to generate the actual number
spike_no = 20;
snr = 40;
% subsamp = 3;
% maxcicles = 100;

% generate the activity dataset with the data given by the bottom topcamera
% calubration and the beads, from the labeled matrix "projectedbeads"
pat = zeros(size(beads_pos, 1), frames);
for idx = 1:size(beads_pos, 1)
%     pat{idx,1} = find(projectedbeads==idx);
    pat(idx,:) = calcium_trace(frames, time_stepsize, spike_no, snr, 'simple');
end
%%
% load the activity from Chen et al, Nature 2013 datasets
% L2/3
% '7' with gcamp6f
% '8' with gcamp6s
% 

% additional spike common in some of the the neurons, to induce crosscorrelation
% across them
% choose a sub set of m neurons, and add to them some spikes in common,
% choice from a poissonian distribution the number of spikes, and choose 5
% neurons which will share these spikes. now start with 3 spikes 

overallcorrelation = 0;
% parameters
frames = 500;
% frames = 2200;

radius_px = 100;
% min_dist_px = 10;
disk_diam_px = 20;
time_stepsize = 0.1;
offset = 1200;

addpath('D:\Users\Comedia\moro\git_code\SpeckledFluoImagingControl\brick-master')
addpath('D:\Users\Comedia\moro\git_code\SpeckledFluoImagingControl\spikes-master')

spikeffold = 'D:\Users\Comedia\moro\git_code\SpeckledFluoImagingControl\spikefinder.train\';
dataset = '8';
% calcium_train = csvread([spikeffold dataset '.train.calcium.csv']);
spike_train = csvread([spikeffold dataset '.train.spikes.csv']);

% generate correlation

for iii = 1:overallcorrelation
    corrneuronsno = 5;
    corrneurons = randsample(size(beads_pos, 1), corrneuronsno);
    interspikeint = 2; % 2 bins between a spike and the other
    spikesno = 3;
    moment = randsample(frames*10, 1); % il 10 e' per il sampling dell'elettrofis
    for nn = corrneurons
        % delete some spikes
        spike_train(randsample(find(spike_train(offset:offset+frames*10,nn)), spikesno), nn) = 0;
        % place them back correlated
        for spp = 1:spikesno
            spike_train(offset+moment+spp*interspikeint,nn) = 1;
        end
    end
end

pat = zeros(size(beads_pos, 1), frames);

% to enhance the crosscorrelation, define a priori a series of poissonian
% events which are in common across all the neurons
% in all the spike_train entries place some common spikes
% for i = <somespikes>
%     spike_train(i, :) = 1
%     %pseudocode


for idx = 1:size(beads_pos, 1)
    % zscore(calcium_train(:,n)) % 100Hz sampling rate, each bin 10 ms
    % find the right scaling of the original dataset, raw, normalized,
    % zscore, or so on...
    % decimate from 100Hz to 10Hz
%     
    % generate calcium activity using the code from Vanzetta group of
    % MLSpike
    spikes = spike_train(offset:frames*10+offset-1,idx);
    % check if there are spikes
    [F F0 drift] = spk_calcium(spikes, 0.01);
    catrace = F(1:frames*10);
    catrace = (catrace - min(catrace(:)))/(max(catrace(:)-min(catrace(:))));
    catrace = catrace(1:10:end);
    pat(idx,:) = catrace;
end
set_bg = 0.0;
%% generate the BG activity 
set_bg = 0.0;
bg_trace = sum(pat,1);
bg_trace = (bg_trace/max(bg_trace(:)))*set_bg;
%% plot the traces to see what happened
figure(); hold;
for idx = 1:20
    plot(pat(idx,:)+idx);
end

if set_bg>0
    plot(bg_trace+21);
end
%% superimpose the whole activity and allow user to select the ROI
%
% topcamera.setDefaultROI();
% pause(2);
% image = zeros(DMD_SIZE);
% for spot = 1:diskssno
%     image(pat{spot,1}) = 255;
% end

% PCOexp = 200;
PCOexp = 400;

PCOsrc.E2ExposureTime = PCOexp;

fullscreen(drawpoints(DMD_SIZE, beads_pos, disk_diam_px), 2);
%
data = getsnapshot(PCOvid);
pause(PCOexp/1000);
fullscreen(rot90(zeros(DMD_SIZE)),2);

figure(10);
% colormap('hsv');
imagesc(data);
h = imrect;
crop = round(h.getPosition);
frame_size = [crop(3),crop(4)];
crop(3:4) = crop(3:4)+crop(1:2);
imagesc(data(crop(2):crop(4),crop(1):crop(3)));

%% acquire data
%
tic
PCOexp = 100;
% PCOexp = 100;

PCOsrc.E2ExposureTime = PCOexp;

% before acquiring data the neurons are singularly projected, to have
% a groundtruth for the spatial decomposition, in order to evaluate the
% decomposition algorithms

% search for a free seed for the dataset, from 1 to 1000 max
for ii = 1:1000
    tmpst = strcat(root_folder,'data_',datestr(now,'ddmmyyyy'),'_',sprintf('%03d',ii) ,'.mat');
    if (exist(tmpst, 'file')~=2)
%         tiff_file = tmpst;
        video_mat_file = tmpst;
        mat_file = strcat(root_folder,'data_',datestr(now,'ddmmyyyy'),'_',sprintf('%03d',ii) ,'_gt.mat');
        break;
    end
end
if ii==1000
    error('please, empty %s',root_folder);
end

main_fig = figure('name','topcamera::dmd');
% better manage the fit of the image
% topcamera_frame = imagesc(rot90(zeros(ceil(frame_size/subsamp))));
topcamera_frame = imagesc(rot90(zeros(ceil(frame_size))));
text1 = text(20,40,'','Color','white');
topcamera_frame.CDataMapping = 'scaled';
video_data = [];

for pt = beads_pos'
    fullscreen(drawpoints(DMD_SIZE, pt, disk_diam_px), 2);
    data = getsnapshot(PCOvid); 
    data = data(crop(2):crop(4),crop(1):crop(3));
    set(topcamera_frame,'CData',data);
    %set(text1,'String',sprintf('grounftruth %d/%d\n',spot,size(beads_pos, 1)));
%     imwrite(data, tiff_file, 'writemode', 'append')
    video_data = cat(3, video_data, data);
end

% for spot = 1:diskssno
%     image = projectedbeads==spot;
%     fullscreen(image,2);
%     data = topcamera.acquireOne();
%     pause(0.1);
%     fullscreen(rot90(zeros(DMD_SIZE)),2);
%     
%     data = data(crop(2):crop(4),crop(1):crop(3));
%     set(topcamera_frame,'CData',data);
%     set(text1,'String',sprintf('grounftruth %d/%d\n',spot,size(beads_pos, 1)));
%     imwrite(data, tiff_file, 'writemode', 'append')
% end


for frame = 1:frames
    image = rot90(zeros(DMD_SIZE));
    for spot = 1:size(beads_pos, 1)
        image = image+ drawpoints(DMD_SIZE, beads_pos(spot, :), disk_diam_px)*pat(spot,frame);
%         image(pat{spot,1}) = image(pat{spot,1}) + pat{spot,2}(frame);
    end
    if set_bg>0
        image = image + drawpoints(DMD_SIZE, bg_beads_pos, disk_diam_px)*bg_trace(frame);
    end
    fullscreen(image,2);
%   facendo cosi' non e' che si sovrappongano continuamente le immagini
%   senza sovrascrivere le vecchie? riempiendo la memoria e cosi' via'...?
    data = getsnapshot(PCOvid);
    pause(0.1);
    fullscreen(rot90(zeros(DMD_SIZE)),2);    
%     data = data(xmin:xmax,ymin:ymax);
    data = data(crop(2):crop(4),crop(1):crop(3));
%     data = data(min(col):max(col),min(row):max(row));
%     imshow(rot90(data(1:3:end,1:3:end)));
%     temp = pre.apply_donut(data);
%     temp = pre.auto_correlation(temp);
    set(topcamera_frame,'CData',data);
%     set(topcamera_frame,'CData',temp(1:subsamp:end,1:subsamp:end));
    set(text1,'String',sprintf('frame %d/%d\n',frame,frames));
%     writeVideo(video,uint8(data/256));
%     imwrite(data, tiff_file, 'writemode', 'append')
    video_data = cat(3, video_data, data);
    % normalizzare meglio, che paiono thresholdate...comunque almeno
    % funziona...
%     refresh;
    drawnow;
    % maybe is not needed the pause for the frame acquisition
%     pause(exposure/1000);
    pause(0.100);
end

% close(video);
fullscreen(rot90(zeros(DMD_SIZE)),2);


save(mat_file,'pat');
save(video_mat_file,'video_data','-v7.3');
% save(video_mat_file,'video_data');
% dmd1.display(rot90(zeros(DMD_SIZE)));
disp('done');
%%
% close everything
delete(dmd1);
imaqreset;
% delete(stage_x);
% delete(stage_y);
closescreen();
clear all