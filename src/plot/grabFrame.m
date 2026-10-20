function img = grabFrame(fig)
%GRABFRAME RGB image of a (possibly invisible) figure, at reduced
%   resolution to keep rendering fast.
img = print(fig, '-RGBImage', '-r72');  % 1280x720 figure -> 960x540 frames
end
