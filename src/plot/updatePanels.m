function updatePanels(pan, t, T, label, ~, ~)
%UPDATEPANELS Per-frame update of the shared video panels.
tr = round(t, 1);                      % avoid 59.96 s printing as 60.0
hh = floor(tr/3600); mm = floor(mod(tr, 3600)/60); ss = mod(tr, 60);
pan.time.String = sprintf('t = %7.1f s (%02d:%02d:%04.1f)   %s', t, hh, mm, ss, label);
pan.dot.XData = t;
pan.dot.YData = interp1(pan.tHist, pan.thetaDeg, t);
pan.curR.Value = t;
pan.curT.Value = t;
setSpacecraftPose(pan.sc, eye(3), zeros(3, 1), T, 3);
end
