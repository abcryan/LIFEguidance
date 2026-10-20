function R = rot2(a)
%ROT2 Active rotation about e2 by angle a, exp(a*[e2]x). Double or CasADi.
c = cos(a); s = sin(a);
R = [c 0 s; 0 1 0; -s 0 c];
end
