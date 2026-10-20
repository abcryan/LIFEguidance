function setSpacecraftPose(h, R, r, T, plumeLen)
%SETSPACECRAFTPOSE Place the spacecraft (rotation R body->scene, position
%   r) and show plumes opposite to each force direction d_k, with length
%   and opacity proportional to T_k/Tmax (plumeLen = length at Tmax, body
%   units).
h.body.Matrix = [h.mag*R, r(:); 0 0 0 1];
for k = 1:numel(T)
    u = T(k)/h.Tmax;
    if u < 1e-3
        h.plume(k).Visible = 'off';
        continue
    end
    ex = -h.thrDir(:,k);                 % exhaust direction
    h.plume(k).Matrix = [alignZ(ex), h.thrPos(:,k); 0 0 0 1] ...
                        * makehgtform('scale', [1 1 max(plumeLen*u, 0.05)]);
    set(h.plume(k).Children, 'FaceAlpha', 0.25 + 0.7*u);
    h.plume(k).Visible = 'on';
end
end

function R = alignZ(v)
z = [0; 0; 1]; v = v/norm(v); c = z'*v;
if c < -1 + 1e-9, R = diag([1 -1 -1]); return, end
K = skew(cross(z, v));
R = eye(3) + K + K*K/(1 + c);
end
