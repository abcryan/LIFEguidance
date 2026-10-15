function log = runClosedLoop(p, ctrl)
%RUNCLOSEDLOOP Apply the NMPC policy (23) to the inertial truth model.
%   The input is recomputed every Ts and held constant in between.
Ts = p.mpc.Ts;
nSteps = round(p.sim.Tend/Ts);
tSub = (0:p.sim.dtLog:Ts)';
nSub = numel(tSub) - 1;
odeOpts = odeset('RelTol', 1e-12, 'AbsTol', 1e-9);

s = initialState(p);
nLog = nSteps*nSub + 1;
log.t = (0:nLog-1)'*p.sim.dtLog;
log.s = zeros(nLog, numel(s));
log.z = zeros(nLog, 39);
log.u = zeros(nLog, 24);
log.mpc = struct('t', {}, 'status', {}, 'iter', {}, 'solveTime', {}, ...
                 'cost', {}, 'Z', {}, 'U', {});
guess = [];

for i = 1:nSteps
    t0 = (i-1)*Ts;
    z = inertialToRelative(s, p);
    [u, info] = ctrl.solve(z, guess);
    fprintf('t = %5.0f s | %-22s | iter %3d | %5.2f s | theta %6.2f deg | |rho| %7.3f m\n', ...
            t0, info.status, info.iter, info.solveTime, ...
            rad2deg(eigenaxisAngle(reshape(z(1:9), 3, 3))), norm(z(13:15)));
    log.mpc(i) = struct('t', t0, 'status', info.status, 'iter', info.iter, ...
                        'solveTime', info.solveTime, 'cost', info.cost, ...
                        'Z', info.Z, 'U', info.U);

    % Shifted warm start for the next solve
    guess.U = [info.U(:,2:end) info.U(:,end)];
    guess.S = [info.S(:,2:end) info.S(:,end)];
    guess.Z = [info.Z(:,2:end) full(ctrl.F(info.Z(:,end), info.U(:,end)))];

    [~, S] = ode113(@(t, x) inertialDynamics(x, u, p), t0 + tSub, s, odeOpts);
    rows = (i-1)*nSub + (1:nSub);
    log.s(rows,:) = S(1:nSub,:);
    log.u(rows,:) = repmat(u', nSub, 1);

    s = S(end,:)';
    s(1:9)   = reshape(projectSO3(reshape(s(1:9), 3, 3)), [], 1);
    s(19:27) = reshape(projectSO3(reshape(s(19:27), 3, 3)), [], 1);
    s(37:39) = min(max(s(37:39), -p.mpc.xiMax), p.mpc.xiMax);
end
log.s(end,:) = s';
log.u(end,:) = log.u(end-1,:);
for j = 1:nLog
    log.z(j,:) = inertialToRelative(log.s(j,:)', p)';
end
log = postProcess(log, p);
end

function log = postProcess(log, p)
n = numel(log.t);
log.theta = zeros(n,1);
log.fc = log.u*p.Bf';
log.tauc = log.u*p.Btau';
log.fNorm = zeros(n,8);
log.alpha1 = nan(n,8);
log.alpha3 = nan(n,8);
log.losViolation = zeros(n,1);
log.gimbalViolation = zeros(n,1);
for j = 1:n
    log.theta(j) = eigenaxisAngle(reshape(log.z(j,1:9), 3, 3));
    log.losViolation(j) = max(p.los.A*log.z(j,13:15)' - p.los.b);
    F = reshape(log.u(j,:), 3, 8);
    log.fNorm(j,:) = vecnorm(F);
    log.gimbalViolation(j) = max(p.thr.A*F, [], 'all');
    for k = 1:8
        if log.fNorm(j,k) > 1e-3
            % (13): f_tk = -R_i(a1) R_k(a3) j-hat |f_tk|
            fk = F(:,k)/log.fNorm(j,k);
            log.alpha3(j,k) = asin(fk(1));
            log.alpha1(j,k) = atan2(-fk(3), -fk(2));
        end
    end
end
log.rho = log.z(:,13:15);
log.rhoDot = log.z(:,16:18);
log.wct = log.z(:,10:12);
end
