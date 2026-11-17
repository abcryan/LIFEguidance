classdef FleetOcpTest < matlab.unittest.TestCase
    %FLEETOCPTEST Problem 4 transcription: the joint NLP (buildFleetOcp) is
    %   the five single-spacecraft Problem 3 blocks (buildOcp, parametric,
    %   swapped with setSpacecraft) coupled only through tau and phiTau.
    %   Coarse mesh (N = 10) for speed; no solve.

    properties
        p
        ref
        fleet
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
