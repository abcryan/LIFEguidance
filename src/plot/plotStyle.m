function s = plotStyle()
%PLOTSTYLE Fixed component colours and pod colours used in every figure.
s.comp = [0.000 0.447 0.741; 0.850 0.325 0.098; 0.466 0.674 0.188];  % x y z
s.pod  = [0.20 0.45 0.80; 0.85 0.35 0.10; 0.30 0.65 0.25; 0.60 0.30 0.70];
s.band = [0.85 0.95 0.85];
% Problems 1, 2 (and 3) in every comparison; Sun-angle violation shading
s.prob = [0.90 0.45 0.05; 0.10 0.40 0.80; 0.45 0.20 0.60];
s.viol = [1.00 0.82 0.80];
% Twenty distinct thruster colours (Tableau-20 order, alternating light/dark)
s.thr = [31 119 180; 255 127 14; 44 160 44; 214 39 40; 148 103 189; 140 86 75; ...
         227 119 194; 127 127 127; 188 189 34; 23 190 207; 174 199 232; 255 187 120; ...
         152 223 138; 255 152 150; 197 176 213; 196 156 148; 247 182 210; ...
         199 199 199; 219 219 141; 158 218 229]/255;
s.figPos = [50 50 1300 850];
end
