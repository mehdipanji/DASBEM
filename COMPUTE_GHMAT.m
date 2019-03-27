                          function  COMPUTE_GHMAT
global X Y CS1 CS2 AT ISTEP NN NE NIE Hd1 Hd2 Gd1 Gd2
%                               FORM  G & H MATRICES
% ASSIGN PRELOCATION FOR H MATRIX
NIN = 2*NIE;
Hd1 = zeros(NIN,NIN);
Gd1 = zeros(NIN,3*NIE);
%
cst = CS1*AT*ISTEP;
%
X1 = X(1:NIN);
Y1 = Y(1:NIN);
%
X1(NIN+1) = X(1);
Y1(NIN+1) = Y(1);
% Calculating the coefficients of H matrix for infinite domain
for ii = 1 : NIN
    for jj = 1 : NIE
        j1 = 2*jj-1;
        j2 = j1+1;
        j3 = j2+1;
        if j3 > NIN
            j3 = 1;
        end
        % AT TIME STEP=1 AND SINGULAR NODES- CALCULATION OF DIAGONAL COEFF.
        if  ((ii-j1)*(ii-j2)*(ii-j3))==0  && ISTEP == 1
            nodo = ii-j1+1;
            if ii == 1 && jj == NIE
                nodo = 3;
            end
            [HT,GT] = SING(X1(ii),Y1(ii),X1(j1),Y1(j1),X1(j2),Y1(j2),X1(j3),Y1(j3),...,
                AT,CS1,nodo);
        else
            % FOR H MATRICES CALCULATE NONSINGULAR COEFF.
            % DIST : DISTANCE FROM COLLOCATOIN NODE TO INTEGRATION ELEMENT
            dist = DISTP(X1(ii),Y1(ii),X1(j1),Y1(j1),X1(j2),Y1(j2),X1(j3),Y1(j3));
            if dist > cst
                continue
            else
                [HT,GT] = NONSING(X1(ii),Y1(ii),X1(j1),Y1(j1),X1(j2),Y1(j2),X1(j3),...,
                    Y1(j3),ISTEP,AT,CS1);
            end
        end
        % ARRANGE G & H COEFF. IN GENERAL G AND H MATRICES
        Hd1(ii,j1) = Hd1(ii,j1) + HT(1);
        Hd1(ii,j2) = Hd1(ii,j2) + HT(2);
        Hd1(ii,j3) = Hd1(ii,j3) + HT(3);
        Gd1(ii,3*jj-2) = GT(1);
        Gd1(ii,3*jj-1) = GT(2);
        Gd1(ii,3*jj) = GT(3);
    end
    if ISTEP == 1
        CI = 0.5;
        if  rem(ii,2) ~= 0 
            i2 = ii-1;
            if ii == 1
                i2 = NIN;
            end
            i1 = i2-1;
            [CI] = HDIAG(X1(i1),Y1(i1),X1(i2),Y1(i2),X1(ii),Y1(ii),X1(ii+1),Y1(ii+1),...,
                X1(ii+2),Y1(ii+2));
        end
        Hd1(ii,ii) = Hd1(ii,ii) + CI;
    end
end
cst = CS2*AT*ISTEP;
%
Hd2 = zeros(NN-NIN,NN-NIN);
Gd2 = zeros(NN-NIN,3*(NE-NIE));
%
X2 = X(NIN+1:end);
Y2 = Y(NIN+1:end);
%
X2(end+1) = X(1);
Y2(end+1) = Y(1);
%
% Calculating the coefficients of H matrix for inclusion
for ii = 1 : (NN-NIN)
    for jj = 1 : (NE-NIE)
        j1 = 2*jj-1;
        j2 = j1+1;
        j3 = j2+1;
        if j3 > (NN-NIN)
            j3 = 1;
        end
        % AT TIME STEP=1 AND SINGULAR NODES- CALCULATION OF DIAGONAL COEFF.
        if  ((ii-j1)*(ii-j2)*(ii-j3))==0  && ISTEP == 1
            nodo = ii-j1+1;
            if ii == 1 && jj == NE-NIE
                nodo = 3;
            end
            [HT,GT] = SING(X2(ii),Y2(ii),X2(j1),Y2(j1),X2(j2),Y2(j2),X2(j3),Y2(j3),...,
                AT,CS2,nodo);
        else
            % FOR H MATRICES CALCULATE NONSINGULAR COEFF.
            % DIST : DISTANCE FROM COLLOCATOIN NODE TO INTEGRATION ELEMENT
            dist = DISTP(X2(ii),Y2(ii),X2(j1),Y2(j1),X2(j2),Y2(j2),X2(j3),Y2(j3));
            if dist > cst
                continue
            else
                [HT,GT] = NONSING(X2(ii),Y2(ii),X2(j1),Y2(j1),X2(j2),Y2(j2),X2(j3),...,
                    Y2(j3),ISTEP,AT,CS2);
            end
        end
        % ARRANGE G & H COEFF. IN GENERAL G AND H MATRICES
        Hd2(ii,j1) = Hd2(ii,j1) + HT(1);
        Hd2(ii,j2) = Hd2(ii,j2) + HT(2);
        Hd2(ii,j3) = Hd2(ii,j3) + HT(3);
        Gd2(ii,3*jj-2) = GT(1);
        Gd2(ii,3*jj-1) = GT(2);
        Gd2(ii,3*jj) = GT(3);
    end
    if ISTEP == 1
        CI = 0.5;
        if  rem(ii,2) ~= 0 
            i2 = ii-1;
            if ii == 1
                i2 = NN-NIN;
            end
            i1 = i2-1;
            [CI] = HDIAG(X2(i1),Y2(i1),X2(i2),Y2(i2),X2(ii),Y2(ii),X2(ii+1),Y2(ii+1),...,
                X2(ii+2),Y2(ii+2));
        end
        Hd2(ii,ii) = Hd2(ii,ii) + CI;
    end
end
% Rearrange the colomns of H mat accoring to the assembling format
Hd3 = Hd2;
k = 0;
for i = 2:NIN
    Hd2(:,i) = Hd3(:,NIN-k);
    k = k+1;
end
%
% Rearrange the colomns of G mat accoring to the assembling format
Gd3 = Gd2;
k = 0;
for i = 1:NIE
    Gd2(:,3*i-2) = Gd3(:,3*(NIE-k));
    Gd2(:,3*i-1) = Gd3(:,3*(NIE-k)-1);
    Gd2(:,3*i-0) = Gd3(:,3*(NIE-k)-2);
    k = k+1;
end
%