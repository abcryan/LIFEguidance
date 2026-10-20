classdef ModelTest < matlab.unittest.TestCase
    %MODELTEST Consistency checks of the parameters, frames and dynamics
    %   (Guidance_Formulation_LIFE.pdf, Secs. 3-4). Run with
    %   runtests('tests') from the project root.

    properties
        p
    end

    methods (TestClassSetup)
        function setup(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            addpath(root); addpath(genpath(fullfile(root, 'src')));
            setupCasadi();
            tc.p = lifeParameters();
        end
    end

    methods (Test)
        function constantsMatchDocument(tc)
            % Digits printed in eqs. (43), (49), (50), (57)
            p = tc.p;
            rel = @(a, b) abs(a - b)/abs(b);
            tc.verifyLessThan(rel(p.eta, 3.00348059454e-6), 1e-11);
            tc.verifyLessThan(rel(p.Om, 1.99098666471e-7), 1e-11);
            tc.verifyLessThan(rel(p.thr.rho, 3.7676252468), 1e-10);
            tc.verifyLessThan(rel(p.beta, 7.61990417335e14), 1e-11);
        end

        function thrusterGeometry(tc)
            p = tc.p;
            tc.verifyEqual(rank(p.B), 6);
            tc.verifyLessThan(norm(p.B*ones(20, 1)), 1e-13);
            tc.verifyEqual(vecnorm(p.thr.pos), 3.8*ones(1, 20), 'AbsTol', 1e-12);
            tc.verifyEqual(vecnorm(p.thr.dir), ones(1, 20), 'AbsTol', 1e-15);
            % Corner cant from the inward radial axis, Sec. 4.3
            corner = setdiff(1:20, 1:5:20);
            a = -p.thr.dir(:, 1:5:20);
            cant = acosd(-sum(p.thr.dir(:, corner).*repelem(a, 1, 4)));
            tc.verifyEqual(cant, repmat(acosd(1/sqrt(3)), 1, 16), 'AbsTol', 1e-10);
            % Positive span: every +-e_i wrench is reachable with T >= 0
            opts = optimoptions('linprog', 'Display', 'off');
            for i = 1:6
                for sgn = [-1 1]
                    e = zeros(6, 1); e(i) = sgn;
                    [~, ~, flag] = linprog(ones(20, 1), [], [], p.B, e, zeros(20, 1), [], opts);
                    tc.verifyEqual(flag, 1, sprintf('wrench %+d e_%d not reachable', sgn, i));
                end
            end
        end

        function frameConventions(tc)
            p = tc.p;
            fr = frames(0, p);
            tc.verifyEqual(fr.Rcb*[1; 0; 0], [1; 0; 0], 'AbsTol', 1e-15);
            % z_b on the ecliptic pole in equatorial axes
            tc.verifyEqual(fr.Rcb(:,3), [0; -sin(p.obliquity); cos(p.obliquity)], 'AbsTol', 1e-15);
            tc.verifyEqual(fr.Rcf*[0; 0; 1], [0; 1; 0], 'AbsTol', 1e-15);   % s^c = e2
            for t = [0 3000 14400]
                fr = frames(t, p);
                tc.verifyLessThan(norm(fr.Rbf'*fr.Rbf - eye(3)), 1e-14);
                % w^f_bf and its derivative against central differences
                dt = 1;
                Rp = frames(t + dt, p).Rbf; Rm = frames(t - dt, p).Rbf;
                W = fr.Rbf'*(Rp - Rm)/(2*dt);       % = [w_bf]x
                % central differences: truncation ~(|w| dt)^2 ~ 1e-8 relative
                tc.verifyEqual([W(3,2); W(1,3); W(2,1)], fr.wbf, 'AbsTol', 1e-8*norm(fr.wbf));
                wd = (frames(t + dt, p).wbf - frames(t - dt, p).wbf)/(2*dt);
                tc.verifyEqual(wd, fr.wbfDot, 'AbsTol', 1e-8*norm(fr.wbfDot));
            end
        end

        function haloPeriod(tc)
            p = tc.p;
            opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-15, 'Events', @crossing);
            [~, ~, te] = ode89(@(~, z) cr3bpNormalized(z, p.eta), [0, 4], p.haloIC, opts);
            T = 2*te(find(te > 0.1, 1))/p.Om/86400;   % skip the crossing at t = 0
            tc.verifyEqual(T, p.haloPeriodDays, 'AbsTol', 1e-4);
            function [v, term, dir] = crossing(t, z)
                v = z(2); term = t > 0.1; dir = 0;
            end
        end

        function gravityDifferenceIsExact(tc)
            p = tc.p;
            R = [1.0083*p.dse; 1e8; 1e8] - p.rbe;
            d = [100; -50; 20];
            exact = double(-p.mue*((vpa(R) + d)/norm(vpa(R) + d)^3 - vpa(R)/norm(vpa(R))^3));
            tc.verifyLessThan(norm(gravityDifference(R, d, p.mue) - exact)/norm(exact), 1e-12);
        end

        function translationMatchesInertialPath(tc)
            % (19)-(22) against the independent inertial-frame route
            p = tc.p; rng(3);
            for k = 1:200
                [t, x, T, rbf] = randomSample(p);
                xd = relativeDynamics(t, x, T, rbf, p);
                a = relativeAccelInertialPath(t, x, T, rbf, p);
                tc.verifyLessThan(norm(xd(11:13) - a), 1e-9*norm(a) + 1e-19);
            end
        end

        function casadiEqualsDouble(tc)
            p = tc.p; rng(4);
            t = casadi.SX.sym('t'); x = casadi.SX.sym('x', 13); T = casadi.SX.sym('T', 20);
            r = casadi.SX.sym('r', 3);
            F = casadi.Function('F', {t, x, T, r}, {relativeDynamics(t, x, T, r, p)});
            for k = 1:50
                [ts, xs, Ts, rs] = randomSample(p);
                a = full(F(ts, xs, Ts, rs)); b = relativeDynamics(ts, xs, Ts, rs, p);
                tc.verifyLessThan(norm(a - b), 1e-14*norm(b));
            end
        end

        function terminalAttitudeIsEquilibrium(tc)
            % R^f_i = I, w^f_fi = 0: the spacecraft spins with the formation
            % about its symmetry axis, which needs no torque
            p = tc.p;
            x = [1; zeros(12, 1)];
            xd = relativeDynamics(1000, x, zeros(20, 1), p.haloIC(1:3)*p.dse, p);
            tc.verifyLessThan(norm(xd(1:7)), 1e-20);
        end

        function quaternionMatchesDcmKinematics(tc)
            % Quaternion form of (10) against the DCM ODE over 4 h
            w = @(t) 1e-3*[sin(1e-3*t); cos(7e-4*t); 0.5];
            qf = @(t, q) 0.5*quatMul([0; w(t)], q);
            Rf = @(t, R) reshape(skew(w(t))*reshape(R, 3, 3), 9, 1);
            o = odeset('RelTol', 1e-13, 'AbsTol', 1e-15);
            R0 = rot2(pi/6);
            [~, Q] = ode89(qf, [0 14400], dcmToQuat(R0), o);
            [~, R] = ode89(Rf, [0 14400], R0(:), o);
            dR = quatToDcm(Q(end,:)')'*reshape(R(end,:), 3, 3);
            tc.verifyLessThan(acos(min((trace(dR) - 1)/2, 1)), 1e-6);
            tc.verifyLessThan(norm(dR - eye(3)), 1e-10);
        end

        function absoluteMatchesRelativePropagation(tc)
            % Same random piecewise-constant thrust in the absolute model
            % (SER + inertial attitude) and the relative model; agreement
            % to 1% of the arrival tolerances over 4 h
            p = tc.p; rng(5);
            N = 48;
            sol.N = N; sol.tNode = linspace(0, p.tf, N + 1);
            sol.T = p.Tmax*rand(20, N).*(rand(20, N) < 0.1)*0.05;
            prop = propagateRelative(sol, p, 2);
            ab = propagateAbsolute(sol, p, prop.t);
            d = stateDiscrepancy(ab.x, prop.x);
            g = p.gate.absRel;
            tc.verifyLessThan(d.max.r, g*p.gate.r);
            tc.verifyLessThan(d.max.v, g*p.gate.v);
            tc.verifyLessThan(d.max.theta, g*p.gate.theta);
            tc.verifyLessThan(d.max.w, g*p.gate.w);
        end
    end
end

function [t, x, T, rbf] = randomSample(p)
t = 14400*rand;
q = randn(4, 1); q = q/norm(q);
x = [q; 1e-3*randn(3, 1); 100*randn(3, 1); 0.1*randn(3, 1)];
T = p.Tmax*rand(20, 1);
rbf = p.haloIC(1:3)*p.dse + 1e7*randn(3, 1);
end
