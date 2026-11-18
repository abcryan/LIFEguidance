classdef FleetOcpTest < matlab.unittest.TestCase
    %FLEETOCPTEST Problem 4 transcription: the joint NLP (buildFleetOcp) is
    %   the five single-spacecraft Problem 3 blocks (buildOcp, parametric,
    %   swapped with setSpacecraft) coupled only through tau and phiTau.
    %   Problem 5: the collision rows (67) are the squared pairwise
    %   separations, the Problem 4 rows are unchanged, and the
    %   closest-approach search of the verification is exact.
    %   Coarse mesh (N = 10) for speed; no solve.

    properties
        p
        ref
        fleet
        fleetCol
        single
    end

    methods (TestClassSetup)
        function setup(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            addpath(root); addpath(genpath(fullfile(root, 'src')));
            setupCasadi();
            p = lifeParameters();
            p.ocp.N = 10;
            p.ocp.ipopt.print_level = 0;
            tc.p = p;
            tc.ref = propagateReference(p, p.tf);
            tc.fleet = buildFleetOcp(p, tc.ref);
            tc.fleetCol = buildFleetOcp(p, tc.ref, [], true);
            tc.single = buildOcp(p, tc.ref, 3, true);
        end
    end

    methods (Test)
        function blocksEqualSingleProblems(tc)
            % g and J of the joint NLP at assembled single-spacecraft
            % points equal the stacked single-spacecraft values
            import casadi.*
            fo = tc.fleet; so = tc.single;
            gF = Function('gF', {fo.nlp.x}, {fo.nlp.g, fo.nlp.f});
            gS = Function('gS', {so.nlp.x, so.nlp.p}, {so.nlp.g, so.nlp.f});
            tauH = 0.9; phi = deg2rad(-20);
            w = zeros(fo.nW, 1); gRef = cell(fo.nS, 1); Jref = 0;
            for i = 1:fo.nS
                pk = spacecraftParams(tc.p, i);
                o = setSpacecraft(so, pk);
                s = unpackW(o, initialGuess(o, pk, phi, (-1)^i), pk);
                U = s.U + 0.02*i;            % distinct per spacecraft
                ws = packW(o, tauH, s.X, s.XC, U, phi);
                w = packW(fo.view(i), tauH, s.X, s.XC, U, phi, w);
                [gi, Ji] = gS(ws, o.par);
                gRef{i} = full(gi); Jref = Jref + full(Ji);
            end
            [g, J] = gF(w);
            g = full(g);
            for i = 1:fo.nS
                tc.verifyEqual(g(fo.view(i).iG), gRef{i}, 'AbsTol', 1e-12);
            end
            tc.verifyEqual(full(J), Jref, 'RelTol', 1e-12);
            tc.verifyEqual(numel(g), fo.nS*so.nG);
        end

        function viewsRoundTrip(tc)
            fo = tc.fleet;
            w = randn(fo.nW, 1);
            for i = 1:fo.nS
                pk = spacecraftParams(tc.p, i);
                s = unpackW(fo.view(i), w, pk);
                w2 = packW(fo.view(i), s.tau/tc.p.tf, s.X, s.XC, s.U, s.phiTau, w);
                tc.verifyEqual(w2, w, 'AbsTol', 1e-12);
            end
        end

        function phaseSensitivityIsTerminalOnly(tc)
            % phiTau enters only the terminal rows (39) of every block
            import casadi.*
            fo = tc.fleet;
            JT = jacobian(fo.nlp.g, fo.nlp.x);
            Jg = Function('Jg', {fo.nlp.x}, {JT(:, 2)});
            nz = find(full(Jg(0.1*randn(fo.nW, 1))));
            tc.verifyTrue(all(ismember(nz, fo.iTermAll)));
            % the combiner's translational target is phase independent:
            % only its attitude rows (1-3) depend on phiTau
            v = fo.view(5);
            d = full(fo.dgTdphi(0.1*randn(fo.nW, 1)));
            dComb = d(ismember(fo.iTermAll, v.iTerm));
            tc.verifyEqual(dComb(4:12), zeros(9, 1), 'AbsTol', 1e-14);
        end

        function collisionRowsAreSquaredSeparations(tc)
            % Problem 5 adds 10 pairs x N*d rows |r_i - r_j|^2/dscale^2 of
            % the collocation positions; all other rows are Problem 4's
            import casadi.*
            fo = tc.fleet; fc = tc.fleetCol;
            tc.verifyEmpty(fo.iCol);
            nC = fc.N*fc.d; nP = size(fc.pairs, 1);
            tc.verifyEqual(nP, 10);
            tc.verifyEqual(numel(fc.iCol), nP*nC);
            tc.verifyEqual(fc.nW, fo.nW);
            tc.verifyEqual(fc.nG, fo.nG + nP*nC);
            w = 0.5*randn(fc.nW, 1);
            gC = Function('gC', {fc.nlp.x}, {fc.nlp.g});
            g4F = Function('g4', {fo.nlp.x}, {fo.nlp.g});
            g = full(gC(w)); g4 = full(g4F(w));
            tc.verifyEqual(g(1:fo.nG), g4, 'AbsTol', 1e-12);
            pos = cell(1, fc.nS);
            for i = 1:fc.nS
                XC = reshape(w(fc.view(i).ix.XC), fc.nx, []).*fc.S;
                pos{i} = XC(8:10, :);
            end
            for k = 1:nP
                dd = vecnorm(pos{fc.pairs(k,1)} - pos{fc.pairs(k,2)}).^2'/fc.colScale^2;
                tc.verifyEqual(g(fc.iCol(fc.colPair == k)), dd, 'RelTol', 1e-12);
            end
            % the bound sets |r_i - r_j| >= dmin + margin, nothing else
            f2 = setCollision(fc, 20, 0.5);
            tc.verifyEqual(f2.lbg(fc.iCol), ((20.5)/fc.colScale)^2*ones(nP*nC, 1), 'RelTol', 1e-15);
            tc.verifyEqual(f2.lbg(1:fo.nG), fo.lbg);
            tc.verifyEqual(f2.ubg, fc.ubg);
        end

        function collisionRowsForListedPairsOnly(tc)
            % pairs = "crossing": only the listed pairs, same rows as in
            % the full version
            import casadi.*
            fc = tc.fleetCol;
            fs = buildFleetOcp(tc.p, tc.ref, [], [1 3; 2 4]);
            nC = fs.N*fs.d;
            tc.verifyEqual(fs.pairs, [1 3; 2 4]);
            tc.verifyEqual(numel(fs.iCol), 2*nC);
            tc.verifyEqual(fs.nG, tc.fleet.nG + 2*nC);
            w = 0.5*randn(fs.nW, 1);
            gS = Function('gS', {fs.nlp.x}, {fs.nlp.g}); gC = Function('gC', {fc.nlp.x}, {fc.nlp.g});
            g = full(gS(w)); gAll = full(gC(w));
            for k = 1:2
                kk = find(ismember(fc.pairs, fs.pairs(k,:), 'rows'));
                tc.verifyEqual(g(fs.iCol(fs.colPair == k)), gAll(fc.iCol(fc.colPair == kk)), 'RelTol', 1e-14);
            end
        end

        function closestApproachStraightLines(tc)
            % two constant-velocity tracks: analytic minimum distance/time
            r0 = [30; -12; 5]; v = [-0.04; 0.01; 0.002];
            rv = @(t) [r0 + v*t; v];
            tStar = -(r0'*v)/(v'*v);
            dStar = norm(r0 + v*tStar);
            [tc_, dc] = closestApproach(rv, 0, 2000);
            tc.verifyEqual(tc_, tStar, 'AbsTol', 1e-6);
            tc.verifyEqual(dc, dStar, 'AbsTol', 1e-9);
            % receding at both ends: the minimum is at the start
            [t0, d0] = closestApproach(rv, 1500, 2000);
            tc.verifyEqual(t0, 1500); tc.verifyEqual(d0, norm(r0 + v*1500), 'AbsTol', 1e-12);
        end

        function formationsSatisfyDmin(tc)
            % Sec. 5.1: smallest separation 25 m initially, 80 m at arrival
            chk = checkCollisionFeasibility(tc.p);
            tc.verifyEqual(chk.initialMin, 25, 'AbsTol', 1e-9);
            tc.verifyEqual(chk.targetMin, 80, 'AbsTol', 1e-9);
            tc.verifyTrue(chk.ok);
        end

        function setSpacecraftChangesOnlyParameters(tc)
            so = tc.single;
            pk = spacecraftParams(tc.p, 5);
            o = setSpacecraft(so, pk);
            tc.verifyEqual(o.par, [initialRelativeState(pk)./so.S; pk.rfp]);
            tc.verifyEqual(o.lbx, so.lbx); tc.verifyEqual(o.lbg, so.lbg);
            tc.verifyEqual(o.ubg, so.ubg);
            pBad = pk; pBad.m = 2*pk.m;
            tc.verifyError(@() setSpacecraft(so, pBad), 'LIFE:physicsMismatch');
        end
    end
end
