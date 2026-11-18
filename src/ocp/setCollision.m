function fo = setCollision(fo, dmin, margin)
%SETCOLLISION Collision distance of a Problem 5 NLP (buildFleetOcp with
%   collision rows): the lower bound of every collision row is
%   ((dmin + margin)/dscale)^2, i.e. |r_i - r_j| >= dmin + margin at every
%   collocation point. dmin = 0 (and margin = 0) recovers Problem 4.
%   Only the bounds change, the compiled solver is reused.
if nargin < 3, margin = 0; end
assert(fo.collision, 'The NLP was built without collision rows.');
fo.lbg(fo.iCol) = ((dmin + margin)/fo.colScale)^2;
fo.dmin = dmin; fo.margin = margin;
end
