classdef FormationTest < matlab.unittest.TestCase
    %FORMATIONTEST Virtual-formation construction and the initial condition
    %   derived from the initial formation (p.target, p.init).

    properties
        p
    end

    methods (TestClassSetup)
        function setup(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            addpath(root); addpath(genpath(fullfile(root, 'src')));
            tc.p = lifeParameters();
        end
    end

    methods (Test)
        function defaultTargetIsEq53(tc)
            % sky = e2, phase 0, 24 h spin reproduces Rx(-pi/2) Rz(ws t)
            p = tc.p;
            tc.verifyEqual(skyBasis([0; 1; 0]), rot1(-pi/2), 'AbsTol', 1e-15);
            for t = [0 5000 14400]
                tc.verifyEqual(frames(t, p).Rcf, rot1(-pi/2)*rot3(2*pi/86400*t), 'AbsTol', 1e-14);
            end
        end

        function skyBasisIsRotationOntoSky(tc)
            rng(11);
            for k = 1:50
                s = randn(3, 1); s = s/norm(s);
                R = skyBasis(s);
                tc.verifyLessThan(norm(R'*R - eye(3)), 1e-14);
                tc.verifyEqual(det(R), 1, 'AbsTol', 1e-14);
                tc.verifyEqual(R(:,3), s, 'AbsTol', 1e-15);
            end
            tc.verifyEqual(skyBasis([0; 0; -1])*[0; 0; 1], [0; 0; -1], 'AbsTol', 1e-15);
        end

        function quaternionMatchesAttitude(tc)
            p = tc.p;
            for t = [0 3000 14400]
                [R, q] = formationAttitude(p.target, t);
                tc.verifyEqual(quatToDcm(q), R, 'AbsTol', 1e-14);
                [R, q] = formationAttitude(p.init, t);
                tc.verifyEqual(quatToDcm(q), R, 'AbsTol', 1e-14);
            end
        end

        function pyramidDimensions(tc)
            p = tc.p;
            check(tc, pyramidPoints(p.target), 80, 480, 960);
            check(tc, pyramidPoints(p.init), 25, 150, 300);
            % Target point: a corner of the target rectangle
            tc.verifyEqual(p.rfp, [40; 240; 0]);
            tc.verifyEqual(p.viz.combiner, [0; 0; 960]);
            function check(tc, P, n, w, h)
                C = P(:, 1:4);
                tc.verifyEqual(C(3,:), zeros(1, 4));                 % rectangle in z = 0
                tc.verifyEqual(norm(C(:,1) - C(:,2)), n, 'AbsTol', 1e-12);   % along x
                tc.verifyEqual(norm(C(:,2) - C(:,3)), w, 'AbsTol', 1e-12);   % along y
                tc.verifyEqual(mean(C, 2), zeros(3, 1), 'AbsTol', 1e-12);    % centred
                tc.verifyEqual(P(:,5), [0; 0; h]);
            end
        end

        function initialFormationPointsOpposite(tc)
            p = tc.p;
            tc.verifyEqual(formationAttitude(p.init, 0)*[0; 0; 1], ...
                           -formationAttitude(p.target, 0)*[0; 0; 1], 'AbsTol', 1e-15);
        end

        function initialStateFromInitialFormation(tc)
            % Map x_i(0) to the absolute model and check it against the
            % initial pyramid directly, for every point, for the default
            % and for a randomly oriented, non-default initial formation
            p0 = tc.p; rng(12);
            s = randn(3, 1);
            alt = struct('sky', s/norm(s), 'phase0', 0.7, 'rate', -3e-4, ...
                         'nulling', 30, 'imaging', 120, 'height', 250, 'point', 1);
            for variant = 1:2
                pv = p0;
                if variant == 2, pv.init = alt; pv = deriveFormations(pv); end
                for k = 1:5
                    pv.init.point = k;
                    verifyInitialState(tc, pv);
                end
            end
        end

        function initialStateDefaultNumbers(tc)
            % Default: corner 1 of the initial rectangle, 180 deg about y_f
            p = tc.p;
            x0 = initialRelativeState(p);
            [~, th] = eigenAngle(x0(1:4));
            tc.verifyEqual(th, pi, 'AbsTol', 1e-12);
            tc.verifyEqual(x0(8:10) + p.rfp, [-12.5; 75; 0], 'AbsTol', 1e-12);
            % Both frames spin at 24 h about opposite axes: relative rate -2 ws
            tc.verifyEqual(x0(11:13), cross([0; 0; -2*p.ws], [-12.5; 75; 0]), 'AbsTol', 1e-15);
        end
    end
end

function verifyInitialState(tc, p)
x0 = initialRelativeState(p);
tc.verifyEqual(norm(x0(1:4)), 1, 'AbsTol', 1e-15);
vf0 = p.haloIC(4:6)*p.dse*p.Om;
z = relToAbs(0, x0, zeros(3, 1), vf0, p);        % absolute state, shared origin
Rcb = rot1(p.obliquity);                          % R^c_b(0)
Rc0 = formationAttitude(p.init, 0);
P = pyramidPoints(p.init);
% Attitude aligned with the initial frame, co-rotating with it
tc.verifyLessThan(norm(quatToDcm(z(1:4)) - Rc0), 1e-14);
tc.verifyEqual(z(5:7), [0; 0; p.init.rate], 'AbsTol', 1e-18);
% Position: the chosen point of the initial pyramid
d = z(8:10) - z(14:16);                           % r_bi - r_bf in O^b
tc.verifyEqual(Rcb*d, Rc0*P(:, p.init.point), 'AbsTol', 1e-12);
% Inertial velocity relative to the shared origin: w_cf0 x r
% (recovered as a difference of ~300 m/s orbital velocities: rounding ~3e-14)
vIn = Rcb*((z(11:13) - z(17:19)) + cross([0; 0; p.Om], d));
tc.verifyEqual(vIn, cross(Rc0*[0; 0; p.init.rate], Rc0*P(:, p.init.point)), 'AbsTol', 1e-12);
end
