function P = pyramidPoints(spec)
%PYRAMIDPOINTS The five points of a virtual formation in its own frame
%   (3x5) [m]: four collectors at the corners of the nulling x imaging
%   rectangle in the z = 0 plane, ordered (+,+), (-,+), (-,-), (+,-) in
%   (x, y), and the combiner at the apex, spec.height along +z.
a = spec.nulling/2; b = spec.imaging/2;
P = [a -a -a  a 0;
     b  b -b -b 0;
     0  0  0  0 spec.height];
end
