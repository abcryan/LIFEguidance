function setFormationPose(h, R)
%SETFORMATIONPOSE Draw the formation pyramid rotated by R (formation axes
%   to scene axes); the formation origin is the scene origin.
P = R*h.V;
nc = size(P, 2) - 1;
apex = P(:, end);
B = P(:, [1:nc, 1]);
E = reshape([P(:, 1:nc); repmat(apex, 1, nc); nan(3, nc)], 3, []);
set(h.base, 'XData', B(1,:), 'YData', B(2,:), 'ZData', B(3,:));
set(h.edges, 'XData', E(1,:), 'YData', E(2,:), 'ZData', E(3,:));
set(h.axis, 'XData', [0 apex(1)], 'YData', [0 apex(2)], 'ZData', [0 apex(3)]);
set(h.dots, 'XData', P(1,:), 'YData', P(2,:), 'ZData', P(3,:));
end
