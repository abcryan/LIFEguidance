function pP = pointingParams(p)
%POINTINGPARAMS Parameters for which the relative model (relativeDynamics,
%   frames, sunCosine, relToAbs, absToRel) is stated in the non-spinning
%   pointing frame O* (Sec. 4.1): the target formation without spin and
%   phase (R^c_f = R^c_*, w_c* = 0) and the state relative to its origin
%   (r_fp = 0), so x = [q*i; w*_*i; r*_*i; rdot*_*i], eqs. (34)-(38).
pP = p;
pP.target.rate = 0;
pP.target.phase0 = 0;
pP.ws = 0;
pP.phi0 = 0;
pP.rfp = zeros(3, 1);
end
