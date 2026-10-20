function img = grabFrame(fig)
%GRABFRAME RGB image of a (possibly invisible) figure at its pixel size.
img = print(fig, '-RGBImage', '-r0');
end
