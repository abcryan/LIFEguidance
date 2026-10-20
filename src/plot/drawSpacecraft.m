function h = drawSpacecraft(ax, p, mag)
%DRAWSPACECRAFT Spacecraft body (Fig. 3a), body axes and 20 exhaust plumes
%   as hgtransforms. Position/attitude via setSpacecraftPose; body drawn
%   magnified by mag about its centre of mass.
h.body = hgtransform(ax);
h.mag = mag;
col = [0.75 0.75 0.78];
drum(h.body, 3.3, -1.2, 0.8, col);
drum(h.body, 2.3, 0.8, 3.2, [0.85 0.85 0.88]);
for j = 0:3
    a = [cos(j*pi/2); sin(j*pi/2); 0];
    box(h.body, 3.55*a, 0.5, [0.45 0.45 0.5]);
end
L = 6;
cols = [0.85 0.1 0.1; 0.1 0.6 0.1; 0.1 0.2 0.9];
for i = 1:3
    e = L*((1:3)' == i);
    plot3([0 e(1)], [0 e(2)], [0 e(3)], 'Color', cols(i,:), 'LineWidth', 2, 'Parent', h.body);
end
[Xc, Yc, Zc] = cylinder([0.12 0.7], 12);
h.plume = gobjects(p.nThr, 1);
for k = 1:p.nThr
    h.plume(k) = hgtransform(h.body);
    surf(Xc, Yc, Zc, 'Parent', h.plume(k), 'FaceColor', [1 0.55 0.1], ...
         'EdgeColor', 'none', 'FaceAlpha', 0.7);
    h.plume(k).Visible = 'off';
end
h.thrPos = p.thr.pos; h.thrDir = p.thr.dir; h.Tmax = p.Tmax;
end

function drum(parent, r, z0, z1, col)
[X, Y, Z] = cylinder(r, 40);
Z = z0 + (z1 - z0)*Z;
surf(X, Y, Z, 'Parent', parent, 'FaceColor', col, 'EdgeColor', 'none', 'FaceLighting', 'gouraud');
th = linspace(0, 2*pi, 41);
patch(r*cos(th), r*sin(th), z1*ones(size(th)), col, 'Parent', parent, 'EdgeColor', 'none');
patch(r*cos(th), r*sin(th), z0*ones(size(th)), col, 'Parent', parent, 'EdgeColor', 'none');
end

function box(parent, c, s, col)
[V, F] = deal([0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1] - 0.5, ...
              [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8]);
patch('Vertices', c' + s*V, 'Faces', F, 'FaceColor', col, 'EdgeColor', 'k', 'Parent', parent);
end
