function img = grabFrame(fig)
%GRABFRAME Current image of a figure, 1280x720 for the video figures:
%   drawn on screen with getframe when the figure is visible (fast),
%   rendered off-screen otherwise (-batch mode).
if strcmp(fig.Visible, 'on')
    drawnow;
    img = getframe(fig).cdata;
    pos = getpixelposition(fig);
    if size(img, 2) ~= pos(3)           % HiDPI screens capture at 2x
        img = imresize(img, [pos(4) pos(3)]);
    end
else
    img = print(fig, '-RGBImage', '-r96');
end
end
