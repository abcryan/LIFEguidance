function S = skew(x)
%SKEW Cross-product matrix, S*y = cross(x,y). Works for double and CasADi.
S = [   0   -x(3)  x(2);
      x(3)    0   -x(1);
     -x(2)  x(1)    0  ];
end
