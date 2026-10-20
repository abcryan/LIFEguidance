function x0 = initialRelativeState(p)
%INITIALRELATIVESTATE x_i(0) of eqs. (55)-(56), quaternion branch q_w > 0.
x0 = [dcmToQuat(p.ic.Rfi); p.ic.wfi; p.ic.rpi; p.ic.vpi];
end
