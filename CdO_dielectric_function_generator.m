clear all
close all
w=[1000:1:8000];
% w=1500
CdO=eps_CdO1(w);


a(:,1)=real (CdO);
a(:,2)=imag (CdO);
plot (w, real(CdO))
hold on 
plot (w, imag(CdO))

pentration=10000000./w/2/pi./imag (sqrt (CdO));
function eps = eps_CdO1(w)
    q=1.60217662e-19;  
    eps_CdO_inf=5.1;
    m=9.10938e-31;
    eps_o=8.854e-12;
    C_CdO=1.47;
    mo_CdO=0.1;

    % n1 is the carrier concentration in cm-3
    n1=4e20;  % CC in cm^-3
    n1_eff=n1*10^-20;

    m_e=mo_CdO*(1+(2*C_CdO)*((0.19732697^2)/(mo_CdO*510998.5))*(3*3.14*n1_eff*10^8)^(2/3))^(1/2)   %effective mass
%     m_e=3
    mu1=200;    % mobility in cm^2/(V-s)
%     mu1=1*(150*0.1622)*m_e^(-1)    % mobility in cm^2/(V-s)
    wp1=(1000/(2*3.14))*sqrt((n1*q^2)/(m_e*m*eps_o));    %plasma frequency of layer 1
    wp1_rad=6.28*wp1;
    wp_f=wp1/sqrt(2);
    gamma1=(10000/(2*3.14))*(q/(mu1*m_e*m));    %damping of layer 1 in Hz
    gam_prod_surf=mu1*m_e;
    gamma_wn=33.35641*(1e-12)*gamma1;
    w_enz_1=sqrt((wp1^2/eps_CdO_inf)-gamma1^2);
    w_enz_1_wn=w_enz_1*(2.998e+10)^-1;
    w_p1_wn=wp1*(2.998e+10)^-1;

    eps=drude(w,wp1,eps_CdO_inf,gamma1);
end

function eps = drude(w,wp1,eps_CdO_inf,gamma1)
    l=length(w);
    i=1;
    wd=(2.998e+10)*w;
    while i<=l
        eps_1_real(1,i)=eps_CdO_inf-(wp1^2/(wd(1,i)^2+gamma1^2));
        eps_1_imag(1,i)=(gamma1*wp1^2)/(wd(1,i)*(wd(1,i)^2+gamma1^2));

        eps(1,i)=eps_1_real(1,i)+1i*eps_1_imag(1,i);

        i=i+1;
    end
end
