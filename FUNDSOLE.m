function  [hd,gd] = FUNDSOLE(ra,ISTEP,AT,CS)
% COMPUTE FUNDAMENTAL SOLUTION COMPONENTS                 
gd = 0;
hd = 0;
a0 = CS*ISTEP*AT/ra;
if a0 <= 1
    return
else
    rt1 = CS*ISTEP*AT/ra;
    rt2 = CS*(ISTEP-1)*AT/ra;
    rt3 = CS*(ISTEP-2)*AT/ra;
    if rt1 > 1
        gdc1 = ISTEP*acosh(rt1);
        gdr1 = sqrt(ISTEP^2-(ra/(CS*AT))^2);
        hd1 = gdr1;
        gd1 = gdc1-gdr1;
    else
        gd1 = 0;
        hd1 = 0;
    end
    if rt2 > 1
        gdc2 = (ISTEP-1)*acosh(rt2);
        gdr2 = sqrt((ISTEP-1)^2-(ra/(CS*AT))^2);
        hd2 = gdr2;
        gd2 = gdc2-gdr2;
    else
        gd2 = 0;
        hd2 = 0;
    end
    if rt3 > 1
        gdc3 = (ISTEP-2)*acosh(rt3);
        gdr3 = sqrt((ISTEP-2)^2-(ra/(CS*AT))^2);
        hd3 = gdr3;
        gd3 = gdc3-gdr3;
    else
        gd3 = 0;
        hd3 = 0;
    end
    gd = (gd1-2*gd2+gd3)/(2*pi);
    hd = -(hd1-2*hd2+hd3)/(2*pi*ra);
end