function updatePanels(pan, t, T, label, x, p)
%UPDATEPANELS Per-frame update of the shared video panels.
tr = round(t, 1);                      % avoid 59.96 s printing as 60.0
hh = floor(tr/3600); mm = floor(mod(tr, 3600)/60); ss = mod(tr, 60);
pan.time.String = {sprintf('t = %8.1f s (%02d:%02d:%04.1f)', t, hh, mm, ss), char(label)};
[~, th] = eigenAngle(x(1:4));
pan.info.String = {sprintf('|r_pi|   = %9.4f m', norm(x(8:10))), ...
                   sprintf('|v_pi|   = %9.2e m/s', norm(x(11:13))), ...
                   sprintf('theta    = %9.4f deg', rad2deg(th)), ...
                   sprintf('|w_fi|   = %9.2e 1/s', norm(x(5:7))), ...
                   sprintf('sum T    = %9.3f N', sum(T))};
pan.bar.YData = T'/p.Tmax;
pan.cursor.Value = t;
end
