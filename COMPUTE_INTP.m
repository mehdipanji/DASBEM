function   COMPUTE_INTP
global U XI YI L NIE NSTEP X Y AT CS1 UI SEIS start EQ fp t0 amax teta GE1
% COMPUTE DISPLACEMENT AT INTERNAL POINTS
%
NIN = 2*NIE;
%
if L == 0
    return
end
GI = zeros(L,3*NIE,NSTEP);
HI = zeros(L,NIN,NSTEP);
%
% prepare node coordinates file
%
r = zeros(1,L);
rr = zeros(1,L);
for ii = 1 : L
    r(ii) = abs(-sin(teta)*XI(ii)+cos(teta)*YI(ii));
    rr(ii) = abs(-sin(teta)*XI(ii)-cos(teta)*YI(ii));
end
t=zeros(1,NSTEP);
for jj = 1 : NSTEP
    if jj == 1
        t(jj) = 0;
    else
        t(jj) = t(jj-1)+AT;
    end
    landa = CS1/fp;
    for kk = 1 : L
        if r(kk) > CS1*t(jj)
            eq1 = 0;
        else
            alpha_inc = CS1*(t(jj)-t0)-sin(teta)*XI(kk)+cos(teta)*YI(kk);
            cst = ((pi/landa)*alpha_inc)^2;
            eq1 = -amax*(2*cst-1)*exp(-cst);
        end
        if rr(kk) > CS1*t(jj)
            eq2 = 0;
        else
            alpha_ref = CS1*(t(jj)-t0)-sin(teta)*XI(kk)-cos(teta)*YI(kk);
            cst = ((pi/landa)*alpha_ref)^2;
            eq2 = -amax*(2*cst-1)*exp(-cst);
        end
        EQ(jj,kk)=eq1+eq2;
    end
end
%find(r==max(r))
%figure(10)
%hold on
%plot(t,EQ(:,1),'b-')
%return
disp('');
disp('***********************INTERNAL POINTS BEGAN************************');
disp('');
%
X1 = X(1:NIN);
Y1 = Y(1:NIN);
%
X1(NIN+1) = X(1);
Y1(NIN+1) = Y(1);
%
for  ISTEP = 1 : NSTEP
    tic;
    for ii = 1 : L
        for jj = 1 : NIE
            j1 = 2*jj-1;
            j2 = j1+1;
            j3 = j2+1;
            if j3 > NIN
                j3 = 1;
            end
            % FOR G and H MATRICES CALCULATE COEFF.
            % DIST : DISTANCE FROM COLLOCATOIN NODE TO INTEGRATION ELEMENT
            dist = DISTP(XI(ii),YI(ii),X1(j1),Y1(j1),X1(j2),Y1(j2),X1(j3),Y1(j3));
            cst = CS1*AT*ISTEP;
            if dist > cst
                continue
            else
                [HT,GT] = NONSING(XI(ii),YI(ii),X1(j1),Y1(j1),X1(j2),Y1(j2),X1(j3),...,
                    Y1(j3),ISTEP,AT,CS1);
            end
            % ARRANGE G & H COEFF. IN GENERAL G AND H MATRICES
            HI(ii,j1,ISTEP) = HI(ii,j1,ISTEP) + HT(1);
            HI(ii,j2,ISTEP) = HI(ii,j2,ISTEP) + HT(2);
            HI(ii,j3,ISTEP) = HI(ii,j3,ISTEP) + HT(3);
            GI(ii,3*jj-2,ISTEP) = GT(1);
            GI(ii,3*jj-1,ISTEP) = GT(2);
            GI(ii,3*jj,ISTEP) = GT(3);
        end
    end
    if ISTEP > 1
        FI = zeros(L,1);
        if SEIS == 1
            for imat = 1 : ISTEP-1
                Un = U(1:NIN,imat);
                Tn = U(NIN+1:end,imat);
                iut = ISTEP-imat+1;
                Hn = HI(:,:,iut);
                Gn = GI(:,:,iut);
                FI = FI + (1/GE1)*Gn*Tn - Hn*Un;
            end
            FI = FI + EQ(ISTEP,:)';
        end
    end
    for ii = 1 : L
        if ISTEP == 1
            if SEIS == 1
                UI(ii,ISTEP) = (1/GE1)*GI(ii,:,ISTEP)*U(NIN+1:end,ISTEP)-HI(ii,:,ISTEP)...,
                    *U(1:NIN,ISTEP)+ EQ(ISTEP,ii);
            else
                UI(ii,ISTEP) = (1/GE1)*GI(ii,:,ISTEP)*U(NIN+1:end,ISTEP)-HI(ii,:,ISTEP)*...,
                    U(1:NIN,ISTEP);
            end
        else
            UI(ii,ISTEP) = FI(ii,1);
        end
    end
    file = ['Int.St.',num2str(ISTEP),'.txt'];               
    nod = fopen(file,'wt');
    for ii = 1 : L
        fprintf(nod,' %-10s\n',num2str(UI(ii,ISTEP)));
    end
    fclose(nod);
    disp(['TIME STEP ( ',num2str(ISTEP),' ) : ',...,
        num2str(toc),'  sec']);
    disp(['TOTAL TIME UNTIL CURRENT STEP :',num2str(clock-start),'  sec']);
end
%