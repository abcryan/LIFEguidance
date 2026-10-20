function setSpacecraftPose(h, R, r, T, plumeLen)
%SETSPACECRAFTPOSE Place the spacecraft (rotation R from O^i to the scene,
%   position r) and show a plume from each firing nozzle, opposite to its
%   force direction d_k. Brightness is proportional to T_k/Tmax; length
%   grows from 30% (any firing thruster stays visible) to plumeLen at
%   Tmax (body units).
h.body.Matrix = [h.mag*R, r(:); 0 0 0 1];
for k = 1:numel(T)
    u = min(T(k)/h.Tmax, 1);
    if u < 1e-3
        h.plume(k).Visible = 'off';
        continue
    end
    len = plumeLen*(0.3 + 0.7*u);
    h.plume(k).Matrix = [alignZ(-h.thrDir(:,k)), h.exit(:,k); 0 0 0 1] ...
                        * makehgtform('scale', [1 1 len]);
    set(h.plume(k).Children, 'FaceColor', [0.95 0.4 0.0] + u*[0.05 0.55 0.6], ...
        'FaceAlpha', 1);
    h.plume(k).Visible = 'on';
end
end
