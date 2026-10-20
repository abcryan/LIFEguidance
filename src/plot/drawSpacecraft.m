function h = drawSpacecraft(ax, p, mag)
%DRAWSPACECRAFT Spacecraft of Fig. 3 in its body frame O^i: upper and lower
%   cylinders, the central ring carrying the four pods on +-x^i, +-y^i,
%   all 20 nozzles at r^i_itk pointing along their exhaust direction -d_k,
%   body axes, and one (hidden) plume per thruster. Position and attitude
%   are set with setSpacecraftPose; the body is drawn magnified by mag.
h.body = hgtransform(ax);
h.mag = mag;
cylinderZ(h.body, 2.5, 0.5, 2.6, [0.86 0.86 0.89]);      % upper cylinder
cylinderZ(h.body, 3.25, -0.5, 0.5, [0.72 0.72 0.76]);    % ring with pods
cylinderZ(h.body, 2.5, -2.6, -0.5, [0.86 0.86 0.89]);    % lower cylinder
for j = 0:3
    a = [cos(j*pi/2); sin(j*pi/2); 0];
    b = [-sin(j*pi/2); cos(j*pi/2); 0];
    podBox(h.body, a, b);
end
nozzleLen = 0.22;
h.exit = zeros(3, p.nThr);
for k = 1:p.nThr
    ex = -p.thr.dir(:,k);               % exhaust direction
    nozzle(h.body, p.thr.pos(:,k), ex, nozzleLen);
    h.exit(:,k) = p.thr.pos(:,k) + nozzleLen*ex;
end
set(findobj(h.body, 'Type', 'surface', '-or', 'Type', 'patch'), ...
    'AmbientStrength', 0.6, 'DiffuseStrength', 0.6, 'SpecularStrength', 0.2);
cols = [0.85 0.1 0.1; 0.1 0.6 0.1; 0.1 0.2 0.9];
for i = 1:3
    e = 5.5*((1:3)' == i);
    plot3([0 e(1)], [0 e(2)], [0 e(3)], 'Color', cols(i,:), 'LineWidth', 2, 'Parent', h.body);
end
% Narrow plumes (half-angle ~4 deg) keep neighbouring thrusters distinct
[Xc, Yc, Zc] = cylinder([0.06 0.2], 10);
h.plume = gobjects(p.nThr, 1);
for k = 1:p.nThr
    h.plume(k) = hgtransform(h.body);
    surf(Xc, Yc, Zc, 'Parent', h.plume(k), 'FaceColor', [1 0.6 0.1], ...
         'EdgeColor', 'none', 'FaceLighting', 'none');
    h.plume(k).Visible = 'off';
end
h.thrDir = p.thr.dir; h.Tmax = p.Tmax;
end

function cylinderZ(parent, r, z0, z1, col)
[X, Y, Z] = cylinder(r, 40);
surf(X, Y, z0 + (z1 - z0)*Z, 'Parent', parent, 'FaceColor', col, 'EdgeColor', 'none');
th = linspace(0, 2*pi, 41);
patch(r*cos(th), r*sin(th), z1*ones(size(th)), col, 'Parent', parent, 'EdgeColor', [0.4 0.4 0.4]);
patch(r*cos(th), r*sin(th), z0*ones(size(th)), col, 'Parent', parent, 'EdgeColor', [0.4 0.4 0.4]);
end

function podBox(parent, a, b)
% Pod housing from the ring surface out to just inside the nozzles
r0 = 3.2; r1 = 3.7; w = 0.5; hz = 0.5;
c = [r0 -w -hz; r1 -w -hz; r1 w -hz; r0 w -hz; r0 -w hz; r1 -w hz; r1 w hz; r0 w hz];
V = c(:,1)*a' + c(:,2)*b' + c(:,3)*[0 0 1];
F = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
patch('Vertices', V, 'Faces', F, 'FaceColor', [0.5 0.5 0.55], 'EdgeColor', 'k', 'Parent', parent);
end

function nozzle(parent, pos, ex, len)
% Bell from the mounting point O^k along the exhaust direction
[X, Y, Z] = cylinder([0.07 0.13], 12);
R = alignZ(ex);
P = R*[X(:)'; Y(:)'; len*Z(:)'] + pos;
surf(reshape(P(1,:), size(X)), reshape(P(2,:), size(X)), reshape(P(3,:), size(X)), ...
     'Parent', parent, 'FaceColor', [0.25 0.25 0.28], 'EdgeColor', 'none');
end
