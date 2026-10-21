function img = grabFrame(fig, height)
%GRABFRAME Current image of a figure at the given video height [px]
%   (width from the figure's aspect ratio): drawn on screen with getframe
%   when the figure is visible (fast; HiDPI captures are resampled),
%   rendered off-screen otherwise (-batch mode).
pos = getpixelposition(fig);
width = 2*round(height*pos(3)/pos(4)/2);  % even, as H.264 requires
if strcmp(fig.Visible, 'on')
    drawnow;
    img = getframe(fig).cdata;
else
    img = print(fig, '-RGBImage', sprintf('-r%d', round(96*height/pos(4))));
end
if size(img, 1) ~= height || size(img, 2) ~= width
    img = imresize(img, [height width]);
end
end
