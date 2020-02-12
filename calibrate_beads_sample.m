% calibrate the tranformation between the camera space and the
% DMD space; then find the beads in the field using the bottom camera, and store the
% positions in a variable int he DMD space to be used as "neurons"

%% here must be changed using the BASLERcamera object
% BASLERvid = videoinput('gentl', 1, 'Mono8');
% BASLERsrc = getselectedsource(BASLERvid);
% BASLERvid.FramesPerTrigger = 1;
%%
BASLERsrc.ExposureTime = 1000;
% check right/left - top/bottom mirroring
sq=120;
calibrationp = checkerboard(sq/2, DMD_SIZE(2)/sq, DMD_SIZE(1)/sq)>0.5;
% maybe instead of the checkerboard will be better to do a simusoidal
% pattern

% put in the center some circles
[XX,YY] = ndgrid(1:DMD_SIZE(2),1:DMD_SIZE(1));
calibrationp((XX-DMD_SIZE(2)/2 -2*sq ).^2 + (YY-DMD_SIZE(1)/2 -2*sq ).^2 <= (sq)^2)=0;
calibrationp((XX-DMD_SIZE(2)/2 -2*sq ).^2 + (YY-DMD_SIZE(1)/2 -2*sq ).^2 <= (sq/2)^2)=255;
calibrationp((XX-DMD_SIZE(2)/2 -2*sq ).^2 + (YY-DMD_SIZE(1)/2 +2*sq ).^2 <= (sq)^2)=0;
calibrationp((XX-DMD_SIZE(2)/2 -2*sq ).^2 + (YY-DMD_SIZE(1)/2 +2*sq ).^2 <= (sq/2)^2)=255;
calibrationp((XX-DMD_SIZE(2)/2 +2*sq ).^2 + (YY-DMD_SIZE(1)/2 -2*sq ).^2 <= (sq)^2)=255;
calibrationp((XX-DMD_SIZE(2)/2 +2*sq ).^2 + (YY-DMD_SIZE(1)/2 -2*sq ).^2 <= (sq/2)^2)=0;

fullscreen(calibrationp,2);
pause(0.4);
fig = figure('position', [200, 200, 500, 500]);
while ishandle(fig),
    calibrationaq = getsnapshot(BASLERvid);
    imshow(calibrationaq);
    drawnow
end
%%
fullscreen(rot90(zeros(DMD_SIZE)),2);
% check if colors hav been acquired, and in that casem, average across
% channels
tmp = size(calibrationaq);
if tmp(end)==3
    calibrationaq = mean(calibrationaq, 3);
end

calibrationaq = calibrationaq-imgaussfilt(calibrationaq, 200);
%%
% first check the periodicity of the pattern
fftim = fft2(calibrationaq);
fftim = fftshift(fftim);
fftim = abs(fftim);
% find the dominant frequencies in the two axis
% the difference of the two means an unespected 
% stretch on the image in one axis
% this must be done in 2d to avoid problems with rotations
[~,locs] = findpeaks(mean(fftim,1),'SortStr','descend');
k1 = (abs(locs(2)-locs(3)))/2;
[~,locs] = findpeaks(mean(fftim,2),'SortStr','descend');
k2 = (abs(locs(2)-locs(3)))/2;
board_size = mean(size(calibrationaq)./[k2, k1]);
% build the cal pattern to find the three spots and calibrate the stuff

calpattern = zeros(size(calibrationaq));

[XX,YY] = ndgrid(1:size(calpattern,1),1:size(calpattern,2));

calpattern((XX-size(calpattern,1)/2).^2 + (YY-size(calpattern,2)/2).^2 <= (board_size)^2) = 1;
calpattern((XX-size(calpattern,1)/2).^2 + (YY-size(calpattern,2)/2).^2 <= (board_size/2)^2) = -1;
cc = normxcorr2(calpattern, calibrationaq);
[ypeak, xpeak] = find(cc==max(cc(:)));
ypeak1 = ypeak-size(calpattern,1)/2;
xpeak1 = xpeak-size(calpattern,2)/2;

calpattern = zeros(size(calibrationaq));
calpattern((XX-size(calpattern,1)/2).^2 + (YY-size(calpattern,2)/2).^2 <= (board_size)^2) = -1;
calpattern((XX-size(calpattern,1)/2).^2 + (YY-size(calpattern,2)/2).^2 <= (board_size/2)^2) = 1;
cc = normxcorr2(calpattern, calibrationaq);
[ypeak, xpeak] = find(cc==max(cc(:)));
ypeak2 = ypeak-size(calpattern,1)/2;
xpeak2 = xpeak-size(calpattern,2)/2;
[XX,YY] = ndgrid(1:size(cc,1),1:size(cc,2));
cc((XX-ypeak).^2 + (YY-xpeak).^2 <= (board_size)^2) = 0;
[ypeak, xpeak] = find(cc==max(cc(:)));
ypeak3 = ypeak-size(calpattern,1)/2;
xpeak3 = xpeak-size(calpattern,2)/2;

% dovrebbe esser giusto, ma meglio controllare sempre
v1 = [xpeak2-xpeak1 ypeak2-ypeak1]';
v2 = [xpeak3-xpeak2 ypeak3-ypeak2]';

v11 = [xpeak3-xpeak1 ypeak3-ypeak1]';
v22 = [xpeak2-xpeak3 ypeak2-ypeak3]';

if dot(v1,v2) > dot(v11,v22)
    spots = [xpeak1 ypeak1; xpeak2 ypeak2; xpeak3 ypeak3];
else
    spots = [xpeak1 ypeak1; xpeak3 ypeak3; xpeak2 ypeak2];
end
    
figure, hold, imagesc(calibrationaq), scatter(spots(:,1), spots(:,2), board_size, [1 0 0], 'filled');
text(spots(:,1)+board_size/3,spots(:,2),num2str((1:3)'), 'Color', [1 0 0], 'FontSize', 20);
%% now I need to find the trasform from the vector "spots"
% into the vector in the DMD space
dmdspots = [DMD_SIZE(1)/2+2*sq DMD_SIZE(2)/2-2*sq; DMD_SIZE(1)/2+2*sq DMD_SIZE(2)/2+2*sq; DMD_SIZE(1)/2-2*sq DMD_SIZE(2)/2+2*sq];

TTspots = [spots'; 1 1 1];

vd1 = dmdspots(2,:)-dmdspots(1,:);
vd2 = dmdspots(3,:)-dmdspots(2,:);
vc1 = spots(2,:)-spots(1,:);
vc2 = spots(3,:)-spots(2,:);

mirror = [sign(dot(vd2,vc2)); sign(dot(vd1,vc1))];
scale = [mirror(1)*norm(vd2)/norm(vc2); mirror(2)*norm(vd1)/norm(vc1)];
Ts = [scale(1) 0 0; 0 scale(2) 0; 0 0 1];
TTspots = Ts*TTspots;

trans = mean(dmdspots)-mean(TTspots(1:2,:)');
Tt = [1 0 trans(1); 0 1 trans(2); 0 0 1];
TTspots = Tt*TTspots;

vc1 = TTspots(1:2,2)'-TTspots(1:2,1)';
vc2 = TTspots(1:2,3)'-TTspots(1:2,2)';
theta = acos(dot(vd2, vc2)/(norm(vd2)*(norm(vc2))));
Tr = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
Tshift = [1 0 -DMD_SIZE(1)/2; 0 1 -DMD_SIZE(2)/2; 0 0 1];
Tbackshift = [1 0 DMD_SIZE(1)/2; 0 1 DMD_SIZE(2)/2; 0 0 1];
TTspots = Tbackshift*Tr*Tshift*TTspots;

TT = Tbackshift*Tr*Tshift*Tt*Ts;
translated_pts = TT*[spots'; 1 1 1];
translated_pts = translated_pts(1:2,:)';

figure, hold, imagesc(calibrationp);
scatter(dmdspots(:,1), dmdspots(:,2), sq*1.2, [0 1 0], 'filled');
text(dmdspots(:,1)-sq,dmdspots(:,2),num2str((1:3)'), 'Color', [0 1 0], 'FontSize', 20);
scatter(translated_pts(:,1), translated_pts(:,2), sq, [1 0 0], 'filled');
text(translated_pts(:,1)+sq/3,translated_pts(:,2),num2str((1:3)'), 'Color', [1 0 0], 'FontSize', 20);


% %% illuminate the whole field and show preview ( FROM NOW USELESS, MOVED TO "theFancyGame...")
% fullscreen(rot90(ones(DMD_SIZE)),2);
% preview(BASLERvid);
% %%
% fullscreen(rot90(ones(DMD_SIZE)),2);
% beads = getsnapshot(BASLERvid);
% fullscreen(rot90(zeros(DMD_SIZE)),2);
% beads = imbinarize(beads, 0.5);
% beads = bwlabel(beads, 8);
% %%
% 
% disks_d = 10;
% projectpattern = zeros(DMD_SIZE(2), DMD_SIZE(1));
% [XX,YY] = ndgrid(1:DMD_SIZE(2),1:DMD_SIZE(1));
% 
% for idx = 1:max(beads(:))
%     [y, x] = find(beads==idx);
%     x = mean(x);
%     y = mean(y);
%     v = TT*[x y 1]';
%     projectpattern((XX-v(2)).^2 + (YY-v(1)).^2 <= (disks_d)^2) = idx;
% end
% fullscreen(projectpattern, 2);
% data = getsnapshot(BASLERvid);
% fullscreen(rot90(zeros(DMD_SIZE)),2);
% 
% figure; imagesc(projectpattern);
% figure; imagesc(beads);
% figure; imagesc(data);
% 
% %% %%%%%%%%%%%%%%%%%%
% % the fancy game of sreaching beads
% 
% disks_d = 10;
% % FOV in percentage of the DMD
% fov = 0.5;
% 
% [XX,YY] = ndgrid(1:DMD_SIZE(1),1:DMD_SIZE(2));
% illumination = (XX-DMD_SIZE(1)/2).^2 + (YY-DMD_SIZE(2)/2).^2 <= (DMD_SIZE(2)*fov*0.5)^2;
% 
% se = strel('disk',5);
% fig = figure('position', [200, 200, 1000, 500]);
% 
% while ishandle(fig),
%     fullscreen(rot90(illumination),2);
%     pause(0.1);
%     orbeads = getsnapshot(BASLERvid);
%     fullscreen(rot90(zeros(DMD_SIZE)),2);
%     beads = imbinarize(orbeads, 0.5);
%     beads = imerode(beads, se);
%     beads = bwlabel(beads, 8);
%     shuffledlbl = randperm(max(beads(:)));
%     tmp = zeros(size(beads));
%     % shuffle the labels
%     for lbl = 1:max(beads(:))
%         tmp(beads==lbl) = shuffledlbl(lbl);
%     end
%     beads = tmp;
%     subplot(1,2,1);
%     imagesc(orbeads);
%     pbaspect([1 1 1]);
%     subplot(1,2,2);
%     hold on
%     imagesc(beads);
%     pbaspect([1 1 1]);
%     for idx = 1:max(beads(:))
%         [y, x] = find(beads==idx);
%         x = mean(x);
%         y = mean(y);
%         text(x, y, int2str(idx), 'Color', [1 0 0], 'FontSize', 15);
%     end
%     hold off
%     
%     drawnow
% end
% 
% projectedbeads = zeros(DMD_SIZE(2), DMD_SIZE(1));
% [XX,YY] = ndgrid(1:DMD_SIZE(2),1:DMD_SIZE(1));
% 
% for idx = 1:max(beads(:))
%     [y, x] = find(beads==idx);
%     x = mean(x);
%     y = mean(y);
%     v = TT*[x y 1]';
%     projectedbeads((XX-v(2)).^2 + (YY-v(1)).^2 <= (disks_d)^2) = idx;
% end
% 
% preview(BASLERvid);
% for iii = 1:3
%     for idx = 1:max(projectedbeads(:))
%         fullscreen((projectedbeads==idx)*255, 2);
%         pause(0.1);
%     end
% end
% fullscreen(rot90(zeros(DMD_SIZE)),2);