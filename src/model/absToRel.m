function x = absToRel(t, z, p)
%ABSTOREL Relative state (33) from the absolute state of absoluteDynamics.
fr = frames(t, p);
qcf = formationQuat(t, p);
qfi = quatMul([qcf(1); -qcf(2:4)], z(1:4));   % R^f_i = (R^c_f)' R^c_i
wfi = quatToDcm(qfi)*z(5:7) - fr.wcf;          % (11)
rfiF = fr.Rbf'*(z(8:10) - z(14:16));
vfiF = fr.Rbf'*(z(11:13) - z(17:19)) - cross(fr.wbf, rfiF);
x = [qfi; wfi; rfiF - p.rfp; vfiF];
end
