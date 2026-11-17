function img = grabFrame(fig, height)
%GRABFRAME Current image of a figure at the given video height [px]
%   (width from the figure's aspect ratio): captured on screen with
%   getframe when the figure is visible and the capture is large enough
%   (fast; HiDPI captures are resampled down), otherwise rendered
%   off-screen at the target resolution (-batch mode, or a video taller
%   than the window), so frames are never upscaled.
pos = getpixelposition(fig);
width = 2*round(height*pos(3)/pos(4)/2);  % even, as H.264 requires
img = [];
if strcmp(fig.Visible, 'on')
    drawnow;
    img = getframe(fig).cdata;
    if size(img, 1) < height, img = []; end
end
if isempty(img)
    img = print(fig, '-RGBImage', sprintf('-r%d', round(96*height/pos(4))));
end
if size(img, 1) ~= height || size(img, 2) ~= width
    img = imresize(img, [height width]);
end
end
