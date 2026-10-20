function R = rot1(a)
%ROT1 Active rotation about e1 by angle a, exp(a*[e1]x). Double or CasADi.
c = cos(a); s = sin(a);
R = [1 0 0; 0 c -s; 0 s c];
end
