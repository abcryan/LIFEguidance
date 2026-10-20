function R = rot3(a)
%ROT3 Active rotation about e3 by angle a, exp(a*[e3]x). Double or CasADi.
c = cos(a); s = sin(a);
R = [c -s 0; s c 0; 0 0 1];
end
