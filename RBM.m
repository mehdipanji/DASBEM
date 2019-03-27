function [HST] = RBM(x,y,N)
% RIGID BODY EFFECT - CALCULATION OF DIAGONAL COEFF. OF H MAT
HST = zeros(N);
%
% N = TOTAL NUMBER OF NODES
% ra = RADIUS
% rd1,rd2,rdn = RADIUS DERIVATIVES
% eta1,eta2 = COMPONENTS OF THE UNIT NORMAL TO THE ELEMENT
% xco,yco = INTEGRATION POINT ALONG THE ELEMENT
% ja = JACOBIAN
%
kes(1) = 0.9739065285;       kes(2) = -0.9739065285;
kes(3) = 0.8650633666;       kes(4) = -0.8650633666;
kes(5) = 0.679409568;        kes(6) = -0.679409568;
kes(7) = 0.4333953941;       kes(8) = -0.4333953941; 
kes(9) = 0.1488743389;       kes(10) = -0.1488743389;
%
w(1) = 0.0666713443;         w(2) = 0.0666713443;
w(3) = 0.1494513491;         w(4) = 0.1494513491;
w(5) = 0.2190863625;         w(6) = 0.2190863625;
w(7) = 0.2692667193;         w(8) = 0.2692667193;
w(9) = 0.2955242247;         w(10) = 0.2955242247;
%
for    ii = 1:N
    for   jj = 1:2:N-1
        kk = jj+2;
        if kk > N
            kk = 1;
        end
        xp = x(ii);
        yp = y(ii);
        x1 = x(jj);
        y1 = y(jj);
        x2 = x(jj+1);
        y2 = y(jj+1);
        x3 = x(kk);
        y3 = y(kk);
        %
        a = x3 - 2*x2 + x1;
        b = (x3-x1)/2;
        c = y3 - 2*y2 + y1;
        d = (y3-y1)/2;
        %
        HW = zeros(1,3);
        %
        if    ii == jj  ||  ii == jj+1 || ii == kk
            e11 = -1;
            e22 = 1;
            nsub = 100;
            deltae = (e22-e11)/100;
            ee1 = e11;
            for k = 1:nsub
                ee2 = ee1 + deltae;
                ja2 = 0.5*(ee2-ee1);
                for i = 1 : 10
                    ee = 0.5*(1-kes(i))*ee1+0.5*(1+kes(i))*ee2;
                    f1 = ee*(ee-1)*0.5;
                    f2 = 1-ee^2;
                    f3 = ee*(ee+1)*0.5;
                    xco = (x1)*f1+(x2)*f2+(x3)*f3;
                    yco = (y1)*f1+(y2)*f2+(y3)*f3;
                    ja1 = sqrt((ee*a+b)^2+(ee*c+d)^2);
                    eta1 = (c*ee+d)/ja1;
                    eta2 = -(a*ee+b)/ja1;
                    % real part
                    ra = sqrt((xco-xp)^2+(yco-yp)^2);
                    rd1 = (xco-xp)/ra;
                    rd2 = (yco-yp)/ra;
                    rdn = rd1*eta1+rd2*eta2;
                    coef = ja1*ja2*w(i);
                    hstr = 1/(2*pi*ra);
                    HW(1,1) = HW(1,1) - hstr*coef*rdn*f1;
                    HW(1,2) = HW(1,2) - hstr*coef*rdn*f2;
                    HW(1,3) = HW(1,3) - hstr*coef*rdn*f3;
                end
                ee1 = ee2;
            end
        else
            for  i = 1:10
                f1 = kes(i)*(kes(i)-1)*0.5;
                f2 = 1-(kes(i)^2);
                f3 = kes(i)*(kes(i)+1)*0.5;
                % real part
                xco = x1*f1 + x2*f2 + x3*f3;
                yco = y1*f1 + y2*f2 + y3*f3;
                ja = sqrt((kes(i)*a+b)^2+(kes(i)*c+d)^2);
                eta1 = (kes(i)*c+d)/ja;
                eta2 = -(kes(i)*a+b)/ja;
                ra = sqrt((xp-xco)^2+(yp-yco)^2);
                rd1 = (xco-xp)/ra;
                rd2 = (yco-yp)/ra;
                rdn = rd1*eta1+rd2*eta2;
                cof = (1/(2*pi))*rdn*w(i)*ja/ra;
                %    
                HW(1,1) = HW(1,1) - cof*f1;
                HW(1,2) = HW(1,2) - cof*f2;
                HW(1,3) = HW(1,3) - cof*f3;
            end
        end
        % PLUG THE HW MATRICES INTO THE GENERAL HST MATRICE.
        HST(ii,jj) = HST(ii,jj) + HW(1);
        HST(ii,jj+1) = HST(ii,jj+1) + HW(2);
        HST(ii,kk) = HST(ii,kk) + HW(3);
    end
end
%
% COMPUTE THE DIAGONAL COEFICIENTS OF THE H MATRIX.
for    ii = 1:N
    for  jj = 1:N
        if  ii ~= jj
           HST(ii,ii) = HST(ii,ii)-HST(ii,jj);
        end
    end
    if   HST(ii,ii) < 0
        HST(ii,ii) = 1 + HST(ii,ii);
    end
end
%