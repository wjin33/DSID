
function [r_p,r_q,r_eps,r_epsel,r_epsed,r_epsE,r_epsid,r_Omega,r_f,r_sigmaT,r_Yd1,r_energy] = DSID_direct_iteration()
%%          Damage Consitutive Model (one element test)

clc
close all
clear all

global nu0 E0 a1 a2 a3 a4 C0 C1 alpha
global FTOL

%% Input paramaeters for DSID model

E0 = 6.8E10;  % unit=Pa
nu0 = 0.21;

a1 = 1.26E-13;
a2 = 3.94E-11;
a3 = -1.26E-12;
a4 = 2.51E-13;
C0 = 1.1e5;
C1 = 2.2E6;
alpha =0.231;



%% Parameters for  computation
theta=0.5;
FTOL = 1e-6;
Iter=150;

%% Store the parameters
Pars=zeros(1,9);
Pars(1)=a1;
Pars(2)=a2;
Pars(3)=a3;
Pars(4)=a4;
Pars(5)=C0;
Pars(6)=C1;
Pars(7)=alpha;
Pars(8)=E0;
Pars(9)=nu0;

%% Load path illustration:
%  eta = % 0 = iso, 2 = pure shear, 3 = triax , 4 = normal and shear stress

% %% Load path example:
% %  2 pure shear
% eta = [2];
% steps = 1;
% chargEps1 = [0]; % confining strain and deviatoric strain
% chargs = [100e6];   % confining stress and deviatoric stress
% ninc = [1000];
% 

% %% Load path example:
% %  3.1 triaxial compression (given deviatoric strain) following 
% %      a isotropic confining compression
% eta = [0 3];
% steps = 2;
% chargEps1 = [0 0.0025]; % confining strain and deviatoric strain
% chargs = [10e6 0];   % confining stress and deviatoric stress
% ninc = [10 1000];
% 

% %% Load path example:
% %  3.2 triaxial compression (given deviatoric stress) following 
% %      a isotropic confining compression
% eta = [0 3];
% steps = 2;
% chargEps1 = [0 0]; % confining strain and deviatoric strain
% chargs = [10e6 150e6];   % confining stress and deviatoric stress
% ninc = [10 1000];
% 
% 

%% Loading path:
% 4 Full stress matrix is given as the laod
eta = [4 4];
steps = 2;
chargEps1 = [0 0];
chargs = [100e6     0;    % s11 vertical
           10e6     0;    % s22
           10e6     0;    % s33
             0   60e6;    % s12
             0      0;    % s23
             0      0];   % s13
ninc = [1000 1000];



%%         INITIALIZATIONS:
sigmaT = zeros(3,3); % non-load
sigmaT_v = mat2_mat1(sigmaT);
Omega = zeros(3,3); % undamaged
Omega_v = mat2_mat1(Omega);  
[De0,S0] = matDO1(Omega);   % undamaded stiffness

%========================================================================
%         Strain decomposition:
%         epsT = epsel + epsed + epsid
%========================================================================
epsT = zeros(3,3);
epsel = zeros(3,3); 
epsed = zeros(3,3);
epsE = zeros(3,3);
epsid = zeros(3,3);
Yd1 = zeros(3,3);
energy = zeros(1,7);
% energy = [sigma:depsilon sigma:depsE sigma:depsel sigma:depsid Y:dOmega C1:trOmega]
% transfer matrix to vector
epsT_v = mat2_mat1(epsT); 
epsel_v = mat2_mat1(epsel); 
epsed_v = mat2_mat1(epsed);
epsE_v = mat2_mat1(epsE);
epsid_v = mat2_mat1(epsid);
Yd1_v = mat2_mat1(Yd1);

%% Storage
r_p(1) = trace(sigmaT)/3;%+pin;
r_q(1) = sigmaT(1,1)-sigmaT(2,2);
r_sigmaT(1,:) = sigmaT_v;
r_epsel(1,:) = epsel_v;
r_epsed(1,:) = epsed_v;
r_epsid(1,:) = epsid_v;
r_epsE(1,:) = epsE_v;
r_Omega(1,:) = Omega_v;
r_Yd1(1,:) = Yd1_v;
r_energy(1,:) = energy;

% yield function
[r_f(1), Yd1] = fdDP(sigmaT,zeros(3,3),Omega,Pars); %trial test

if r_f(1) > 0
    error('non-elastic initial state')
end


%% SIMULATION (for a single stress path component):


tinc = 1;
for icharg = 1:steps 
    disp(['============= load step #',num2str(icharg),' ============='])
    
    for inc = 1:ninc(icharg) % load increments
        disp(['              increments #',num2str(inc),'              '])
        
        [matDOm,SOm] = matDO1(Omega);% matDOm SOm 4th tensor
        
        if eta(icharg) == 0 %iso
            if chargs(icharg)~=0   % stress controlled
                dsig = chargs(icharg)/ninc(icharg)*eye(3);
                deps = Aijkl_Bkl(SOm,dsig);
            else    % strain controlled
                deps = chargEps1(icharg)/ninc(icharg)*eye(3);
                dsig = Aijkl_Bkl(matDOm,deps);
            end
            
        elseif eta(icharg) == 3 %triaxial
            % Elastic trial:
            if chargs(icharg)~=0   % stress controlled
                dsig = chargs(icharg)/ninc(icharg)*[1 0 0;0 0 0;0 0 0];
                deps = Aijkl_Bkl(SOm,dsig);
            else    % strain controlled
                deps = chargEps1(icharg)/ninc(icharg)*[1 0 0;0 0 0;0 0 0];
                matDOm_2=mat4_mat2(matDOm,1);% matDOm_2 matrix
                
                deps(2,2) = -matDOm_2(2,1)*deps(1,1)/(matDOm_2(2,2)+matDOm_2(2,3));
                deps(3,3) = deps(2,2);
                dsig = Aijkl_Bkl(matDOm,deps);
                dsig(2,2)=0;
                dsig(3,3)=0;
             
            end
        elseif eta(icharg) == 2 % shear
            % Elastic trial:
            if chargs(icharg)~=0   % stress controlled
                dsig = chargs(icharg)/ninc(icharg)*[0 1 0;1 0 0;0 0 0];
                deps = Aijkl_Bkl(SOm,dsig);
            else    % strain controlled
                deps = chargEps1(icharg)/ninc(icharg)*[0 1 0;1 0 0;0 0 0];
                matDOm_2=mat4_mat2(matDOm,1);% matDOm_2 matrix
                
                deps(2,2) = -matDOm_2(2,1)*deps(1,1)/(matDOm_2(2,2)+matDOm_2(2,3));
                deps(3,3) = deps(2,2);
                dsig = Aijkl_Bkl(matDOm,deps);
                dsig(2,2)=0;
                dsig(3,3)=0;
             
            end
        elseif eta(icharg) == 4 % shear and normal stress
             % only stress controlled
             dsig = 1/ninc(icharg)*[chargs(1,icharg) chargs(4,icharg) chargs(6,icharg);
                                    chargs(4,icharg) chargs(2,icharg) chargs(5,icharg);
                                    chargs(6,icharg) chargs(5,icharg) chargs(3,icharg)];
             deps = Aijkl_Bkl(SOm,dsig);
        end
               
        [ fd0,Ydm] = fdDP(sigmaT,zeros(3,3),Omega,Pars);
        [ fd,Yd1] = fdDP(sigmaT,dsig,Omega,Pars);
% 
        nbStep = 0;
        if fd <= FTOL % elastic increment

            sigmaT = sigmaT + dsig;
            sigEner = sigmaT-theta*dsig;
            %deps = Aijkl_Bkl(SOm,dsig);
            depsel = Aijkl_Bkl(S0,dsig);
            depsid = zeros(3,3);
            dOmega = zeros(3,3);
            trdOmega = 0;
            depsE = deps;
            depsed = deps - depsel;

            epsE = epsE + depsE;
            epsT = epsT + deps ;
            epsed = epsed + depsed;
            epsid = epsid + depsid;
            epsel = epsel + depsel; 

            pn = trace(sigmaT)/3;
            qn = sigmaT(1,1) - sigmaT(2,2);
            energy(1,1)=Aij_Bij(sigEner,deps);
            energy(1,2)=Aij_Bij(sigEner,depsE);
            energy(1,3)=Aij_Bij(sigEner,depsel);
            energy(1,4)=Aij_Bij(sigEner,depsid);
            energy(1,5)=Aij_Bij(Yd1,dOmega);
            dG_dtrdO = (a1+a3)*(trace(sigEner))^2+(a2+a4)*trace(sigEner*sigEner);
            trdOmega = trace(dOmega);
            energy(1,6)=dG_dtrdO*trdOmega;
            energy(1,7)=0.5*Aij_Bij(sigmaT,epsE);
            
            end_r = length(r_p)+1;
            r_p(end_r) = pn;
            r_q(end_r) = qn;
            % transfer matrix to vector
            sigmaT_v = mat2_mat1(sigmaT);
            epsT_v = mat2_mat1(epsT);
            epsel_v = mat2_mat1(epsel); 
            epsed_v = mat2_mat1(epsed);
            epsE_v = mat2_mat1(epsE);
            epsid_v = mat2_mat1(epsid);
            Omega_v = mat2_mat1(Omega);
            
            r_sigmaT(end_r,:)=sigmaT_v;
            r_eps(end_r,:) = epsT_v;
            r_epsel(end_r,:) = epsel_v;
            r_epsed(end_r,:) = epsed_v;
            r_epsE(end_r,:) = epsE_v;
            r_epsid(end_r,:) = epsid_v;
            r_Omega(end_r,:) = Omega_v;
            [r_f(end_r),Yd1] = fdDP(sigmaT,zeros(3,3),Omega,Pars);
            Yd1_v = mat2_mat1(Yd1);
            r_Yd1(end_r,:) = Yd1_v;
            r_energy(end_r,1:6) = r_energy(end_r-1,1:6)+energy(1,1:6);
            r_energy(end_r,7) = energy(1,7);
        else          
            depsid = zeros(3,3);
            depsE = zeros(3,3);
            if (chargs(icharg)~=0 || eta(icharg) == 4)
                 incinc=1;
                lambdam=0;lambdan=10^-3;
              while (abs(fd)>FTOL) && (incinc<Iter) %  for incinc=1:50%
               [ fdm,epsm,Omegam,depsidm ] = fd_lam2( lambdam, Omega, sigmaT, epsid, Pars,dsig);
               [ fdn,epsn,Omegan,depsidn ] = fd_lam2( lambdan, Omega, sigmaT, epsid, Pars,dsig);
               lambda = lambdan - fdn*(lambdan-lambdam)/(fdn-fdm);
                [fd,eps,OmegaT,depsidT] = fd_lam2( lambda, Omega, sigmaT, epsid, Pars,dsig);          
                dOmega = OmegaT-Omega;                
                sigmaTT = sigmaT + dsig;
                
                [ fd,Yd1] = fdDP(sigmaTT,zeros(3,3),OmegaT,Pars);
                incinc=incinc+1;
                lambdam=lambdan;
                lambdan=lambda;
              end
                deps = eps - epsT;
                depsid = depsidT;
                depsel = Aijkl_Bkl(S0,dsig);
                depsed = deps-depsid-depsel;
                depsE = depsed + depsel;
                epsT = epsT + deps;
                Omega = Omega + dOmega;

 %              depsed = depsE - depsel;
                epsel = epsel + depsel; 
                epsid = epsid + depsid;
                epsE = epsE + depsE;  
                epsed = epsed + depsed;

            else
                  incinc=1;
                lambdam=0;lambdan=10^-3;
              while (abs(fd)>FTOL) && (incinc<Iter) %  for incinc=1:50%
                [fdm,sigm,Omegam,depsidm] = fd_lam( lambdam, Omega, sigmaT, epsT, epsid, Pars,deps);
                [fdn,sign,Omegan,depsidn] = fd_lam( lambdan, Omega, sigmaT, epsT, epsid, Pars,deps);
                lambda = lambdan - fdn*(lambdan-lambdam)/(fdn-fdm);
                [fd,sig,OmegaT,depsidT] = fd_lam( lambda, Omega, sigmaT, epsT, epsid, Pars,deps);
                sigmaTT = sig;          
                dsig = sigmaTT-sigmaT;
                dOmega = OmegaT-Omega;
                [fd,Yd1] = fdDP(sigmaTT,zeros(3,3),OmegaT,Pars);
                incinc=incinc+1;
                lambdam=lambdan;
                lambdan=lambda;
              end
                          
              depsid = depsidT;
              epsT = epsT + deps;
              Omega = OmegaT;
              depsel = Aijkl_Bkl(S0,dsig);
              depsed = deps - depsel-depsid;
              depsE = depsel+depsed;
              
              epsel = epsel + depsel; 
              epsid = epsid + depsid;
              epsE = epsE + depsE;       
              epsed = epsed + depsed;           
            end
            sigmaT = sigmaT + dsig;
            sigEner = sigmaT-theta*dsig;
            pn = trace(sigmaT)/3;
            qn = sigmaT(1,1) - sigmaT(2,2);
            energy(1,1) = Aij_Bij(sigEner,deps);
            energy(1,2) = Aij_Bij(sigEner,depsE);
            energy(1,3) = Aij_Bij(sigEner,depsel);
            energy(1,4) = Aij_Bij(sigEner,depsid);
            
            dG_dtrdO = (a1+a3)*(trace(sigEner))^2+(a2+a4)*trace(sigEner*sigEner);
            trdOmega = trace(dOmega);
            energy(1,6) = dG_dtrdO*trdOmega;
            energy(1,7)=0.5*Aij_Bij(sigmaT,epsE);
            
            end_r = length(r_p)+1;
            r_p(end_r) = pn;
            r_q(end_r) = qn;
            sigmaT_v = mat2_mat1(sigmaT);
            epsT_v = mat2_mat1(epsT);
            epsel_v = mat2_mat1(epsel); 
            epsed_v = mat2_mat1(epsed);
            epsE_v = mat2_mat1(epsE);
            epsid_v = mat2_mat1(epsid);
            Omega_v = mat2_mat1(Omega);
            
            r_sigmaT(end_r,:) = sigmaT_v;
            r_eps(end_r,:) = epsT_v;
            r_epsel(end_r,:) = epsel_v;
            r_epsed(end_r,:) = epsed_v;
            r_epsE(end_r,:) = epsE_v;
            r_epsid(end_r,:) = epsid_v;
            r_Omega(end_r,:) = Omega_v;

            [fd1,Yd1]=fdDP(sigmaT,zeros(3,3),Omega,Pars);
            [r_f(end_r),Yd1] = fdDP(sigmaT,zeros(3,3),Omega,Pars);
            Yd1_v = mat2_mat1(Yd1);
            r_Yd1(end_r,:)=Yd1_v;
            energy(1,5) = Aij_Bij(Yd1,dOmega);
            r_energy(end_r,1:6) = r_energy(end_r-1,1:6)+energy(1,1:6);
            r_energy(end_r,7) = energy(1,7);

        end
        tinc = tinc + 1;% just counts total increments (if several loadings)
    end % increments
    r_ends(icharg) = length(r_p);
     
end % loading parts

%% POST-PROCESSING:

figure('Name','q(eps1)','NumberTitle','off')
plot(r_eps(:,1),r_q/1000000,'Linewidth',3)
xlabel('Axial strain, \epsilon_1','FontSize',20)
ylabel('Deviatoric stress, q (MPa)','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'q_eps1.eps'

figure1=figure('Name','q(eps1,eps3)','NumberTitle','off');
plot(r_eps(:,1),r_q/1000000,'-b',r_eps(:,3),r_q/1000000,'-r','Linewidth',3)
xlabel('Axial strain, \epsilon_1','FontSize',20)
ylabel('Deviatoric stress, q (MPa)','FontSize',20)
legend('\epsilon_1', '\epsilon_3','Location','Best')
grid
set(gca,'FontSize',20)
set(gca,'XDir','reverse')
set(gca,'YDir','reverse')
print -depsc 'q_eps1-3.eps'

figure('Name','q(eps3)','NumberTitle','off')
plot(r_eps(:,3),r_q/1000000,'Linewidth',3)
xlabel('Axial strain, \epsilon_3','FontSize',20)
ylabel('Deviatoric stress, q (MPa)','FontSize',20)

grid
set(gca,'FontSize',20)
print -depsc 'q_eps3.eps'

figure('Name','Omega(q)','NumberTitle','off')
plot(r_q/1000000,r_Omega(:,1),'-b','Linewidth',3)
hold on
plot(r_q/1000000,r_Omega(:,3),'--r','Linewidth',5)
xlabel('Deviatoric stress, q (MPa)','FontSize',20)
ylabel('Damage variable, D','FontSize',20)
legend('D_{11}', 'D_{33}','Location','Best')
grid
set(gca,'FontSize',20)
print -depsc 'Omega_q.eps'

figure('Name','Omega(eps1)','NumberTitle','off')
plot(r_eps(:,1),r_Omega(:,1),'-b','Linewidth',3)
hold on
plot(r_eps(:,1),r_Omega(:,3),'--r','Linewidth',5)
xlabel('Axial strain, \epsilon_1','FontSize',20)
ylabel('Damage variable, D','FontSize',20)
legend('D_{11}', 'D_{33}','Location','Best')
grid
set(gca,'FontSize',20)
print -depsc 'Omega_eps1.eps'

figure('Name','sigma12(eps12)','NumberTitle','off')
plot(r_eps(:,4),r_sigmaT(:,4)/1000000,'Linewidth',3)
xlabel('Axial strain, \epsilon_{12}','FontSize',20)
ylabel('Deviatoric stress, \sigma_{12} (MPa)','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'sigma12_eps12.eps'

figure('Name','Omega3(eps3)','NumberTitle','off')
plot(r_eps(:,3),r_Omega(:,3),'Linewidth',3)
xlabel('Axial strain, \epsilon_3','FontSize',20)
ylabel('Damage variable, \Omega_3 (-)','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'Omega3_eps3.eps'

figure('Name','Omega1(eps1)','NumberTitle','off')
plot(r_eps(:,1),r_Omega(:,1),'Linewidth',3)
xlabel('Axial strain, \epsilon_1','FontSize',20)
ylabel('Damage variable, \Omega_1 (-)','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'Omega1_eps1.eps'

figure('Name','Omega3(eps1)','NumberTitle','off')
plot(r_eps(:,1),r_Omega(:,3),'Linewidth',3)
xlabel('Axial strain, \epsilon_1','FontSize',20)
ylabel('Damage variable, \Omega_3 (-)','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'Omega3_eps1.eps'

figure('Name','Omega2(eps1)','NumberTitle','off')
plot(r_eps(:,1),r_Omega(:,2),'Linewidth',3)
xlabel('Axial strain, \epsilon_1','FontSize',20)
ylabel('Damage variable, \Omega_2 (-)','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'Omega2_eps1.eps'

figure('Name','Omega1(eps12)','NumberTitle','off')
plot(r_eps(:,4),r_Omega(:,1),'Linewidth',3)
xlabel('Axial strain, \epsilon_{12}','FontSize',20)
ylabel('Damage variable, \Omega_1 (-)','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'Omega1_eps12.eps'

figure('Name','Omega12(eps12)','NumberTitle','off')
plot(r_eps(:,4),r_Omega(:,4),'Linewidth',3)
xlabel('Axial strain, \epsilon_{12}','FontSize',20)
ylabel('Damage variable, \Omega_{12} (-)','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'Omega12_eps12.eps'

figure('Name','f','NumberTitle','off')
plot(r_f,'Linewidth',3)
xlabel('number','FontSize',20)
ylabel('Damage function, f','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'f.eps'

figure1 = figure('NumberTitle','off','Name','Energy1');
% Create axes
axes1 = axes('Parent',figure1,'FontSize',18);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
% Create multiple lines using matrix input to plot
plot1 = plot(r_energy(:,1),'Parent',axes1,'LineWidth',3);
plot2 = plot(r_energy(:,7),'Parent',axes1,'Linewidth',3);
plot3 = plot(r_energy(:,4),'Parent',axes1,'Linewidth',3);
plot4 = plot(r_energy(:,5),'Parent',axes1,'Linewidth',3);

set(plot1,'Color',[0 0 1],'DisplayName','External work');
set(plot2,'LineStyle','--','Color',[1 0 0],...
    'DisplayName','Elastic strain energy');
set(plot3,'LineStyle',':','Color',[0 1 0],...
    'DisplayName','Irreveasible strain energy');
set(plot4,'LineStyle','-.','Color',[0 1 1],...
    'DisplayName','Crack debonding');
% Create xlabel
xlabel('Substeps','FontSize',20);
% Create ylabel
ylabel('Energy','FontSize',20);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','Best','FontSize',15);
print -depsc 'energy1.eps'

figure1 = figure('NumberTitle','off','Name','EnergyPercentage');
% Create axes
axes1 = axes('Parent',figure1,'FontSize',18);
% xlim(axes1,[0 1100]);
% ylim(axes1,[0 100]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
% Create multiple lines using matrix input to plot
end_r = length(r_energy);
plot1 = plot(100.*r_energy(2:end_r,1)./r_energy(2:end_r,1),'Parent',axes1,'LineWidth',3);
plot2 = plot(100.*r_energy(2:end_r,7)./r_energy(2:end_r,1),'Parent',axes1,'Linewidth',3);
plot3 = plot(100.*r_energy(2:end_r,4)./r_energy(2:end_r,1),'Parent',axes1,'Linewidth',3);
plot4 = plot(100.*r_energy(2:end_r,5)./r_energy(2:end_r,1),'Parent',axes1,'Linewidth',3);

set(plot1,'Color',[0 0 1],'DisplayName','External work');
set(plot2,'LineStyle','--','Color',[1 0 0],...
    'DisplayName','Elastic strain energy');
set(plot3,'LineStyle',':','Color',[0 1 0],...
    'DisplayName','Irreveasible strain energy');
set(plot4,'LineStyle','-.','Color',[0 1 1],...
    'DisplayName','Crack debonding');
% Create xlabel
xlabel('Substeps','FontSize',20);
% Create ylabel
ylabel('E_i/E_{Total}, %','FontSize',20);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','Best','FontSize',15);
print -depsc 'energypercentage.eps'
end

%%                  Functions

function [ scalar ] = Aij_Bij( A,B )
%   Double contraction: scalar = A_ij*B_ij

scalar = 0;
for i= 1:3
    for j=1:3 
        scalar = scalar + A(i,j)*B(i,j); 
    end
end
   
end

function [ C ] = Aij_Bkl( A,B )

C(1:3,1:3,1:3,1:3) = 0;
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                C(i,j,k,l) = A(i,j)*B(k,l);
            end
        end
    end
end

end

function [ C ] = Aijkl_Bij( A,B )

C = zeros(3,3);
for k = 1:3
    for l = 1:3
        for i = 1:3
            for j = 1:3
              C(k,l) = C(k,l)+A(i,j,k,l)*B(i,j);
            end
        end
    end
end

end

function [ C ] = Aijkl_Bkl( A,B )

C = zeros(3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
              C(i,j) = C(i,j)+A(i,j,k,l)*B(k,l);
            end
        end
    end
end

end

function [ dY_dsig ] = dY_dsigf( sigmaT )
%
global E0 nu0 a1 a2 a3 a4
E = eye(3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                dY_dsig(i,j,k,l)=2*a1*trace(sigmaT)*E(i,j)*E(k,l)+1/2*a2*(E(i,k)*sigmaT(l,j)+...
                    E(i,l)*sigmaT(j,k)+E(j,l)*sigmaT(i,k)+E(j,k)*sigmaT(i,l))+a3*(E(k,l)*sigmaT(i,j)+...
                1/2*trace(sigmaT)*(E(i,k)*E(j,l)+E(i,l)*E(j,k)))+2*a4*sigmaT(k,l)*E(i,j);
            end
        end
    end
end

end

function [ f,Yd1 ] = fdDP( sigmaT,dsig,Omega,Pars)
%       Damage Yield Function (Pseudo-Drucker-Pager)
a1=Pars(1);
a2=Pars(2);
a3=Pars(3);
a4=Pars(4);
C0=Pars(5);
C1=Pars(6);
alpha=Pars(7);
sigmaTT = sigmaT + dsig ;
trsigmaT = trace(sigmaTT);
trOmega = trace(Omega);
Yd1 = a1*(trsigmaT)^2*eye(3)+a2*sigmaTT*sigmaTT+a3*trsigmaT*sigmaTT+a4*trace(sigmaTT*sigmaTT)*eye(3);
P_1=matP_1(sigmaTT);

P_1Y=Aijkl_Bkl(P_1,Yd1);
trY = trace(P_1Y);
S = P_1Y-1/3*trY*eye(3);
f = sqrt(0.5*Aij_Bij(S,S))-alpha*trY-C0-C1*trOmega;
end

function [ invA ] = invmat4( A )
% calculate inverse of 4th-order tensor
E6=eye(6);
A_2=mat4_mat2(A,2);
invA_2=A_2\E6;
invA=mat2_mat4(invA_2,1);
end

function [ vector ] = mat2_mat1( matrix )
%==========================================================================
%
%    MAT2_MAT1 
%    Transfer a 3*3 matrix to 6*1 vextor (default format in ABAQUS)
%
%    sig_11 sig_12 sig_13
%    sig_21 sig_22 sig_23 ==> [sig_11 sig_22 sig_33 sig_12 sig_13 sig_23]^T
%    sig_31 sig_32 sig_33
%
%==========================================================================
vector = zeros(6,1);
for i = 1:3
    vector(i) = matrix(i,i);
    for j = i+1:3
        vector(i+j+1) = matrix(i,j);
    end
end

end

function [ tensor ] = mat2_mat4( matrix,coe )
      tensor(1:3,1:3,1:3,1:3) = 0;
      if coe==1
          coe1=1;
          coe2=1;
      elseif coe==2
          coe1=2;
          coe2=4;
      end
      tensor(1,1,1,1) = matrix(1,1);
      tensor(1,1,2,2) = matrix(1,2);
      tensor(1,1,3,3) = matrix(1,3);
      tensor(1,1,1,2) = matrix(1,4)/coe1;
      tensor(1,1,2,1) = matrix(1,4)/coe1;
      tensor(1,1,2,3) = matrix(1,5)/coe1;
      tensor(1,1,3,2) = matrix(1,5)/coe1;
      tensor(1,1,1,3) = matrix(1,6)/coe1;
      tensor(1,1,3,1) = matrix(1,6)/coe1;

      tensor(2,2,1,1) = matrix(2,1);
      tensor(2,2,2,2) = matrix(2,2);
      tensor(2,2,3,3) = matrix(2,3);
      tensor(2,2,1,2) = matrix(2,4)/coe1;
      tensor(2,2,2,1) = matrix(2,4)/coe1;
      tensor(2,2,2,3) = matrix(2,5)/coe1;
      tensor(2,2,3,2) = matrix(2,5)/coe1;
      tensor(2,2,1,3) = matrix(2,6)/coe1;
      tensor(2,2,3,1) = matrix(2,6)/coe1;

      tensor(3,3,1,1) = matrix(3,1);
      tensor(3,3,2,2) = matrix(3,2);
      tensor(3,3,3,3) = matrix(3,3);
      tensor(3,3,1,2) = matrix(3,4)/coe1;
      tensor(3,3,2,1) = matrix(3,4)/coe1;
      tensor(3,3,2,3) = matrix(3,5)/coe1;
      tensor(3,3,3,2) = matrix(3,5)/coe1;
      tensor(3,3,1,3) = matrix(3,6)/coe1;
      tensor(3,3,3,1) = matrix(3,6)/coe1;

      tensor(1,2,1,1) = matrix(4,1)/coe1;
      tensor(1,2,2,2) = matrix(4,2)/coe1;
      tensor(1,2,3,3) = matrix(4,3)/coe1;
      tensor(1,2,1,2) = matrix(4,4)/coe2;
      tensor(1,2,2,1) = matrix(4,4)/coe2;
      tensor(1,2,2,3) = matrix(4,5)/coe2;
      tensor(1,2,3,2) = matrix(4,5)/coe2;
      tensor(1,2,1,3) = matrix(4,6)/coe2;
      tensor(1,2,3,1) = matrix(4,6)/coe2;

      tensor(2,3,1,1) = matrix(5,1)/coe1;
      tensor(2,3,2,2) = matrix(5,2)/coe1;
      tensor(2,3,3,3) = matrix(5,3)/coe1;
      tensor(2,3,1,2) = matrix(5,4)/coe2;
      tensor(2,3,2,1) = matrix(5,4)/coe2;
      tensor(2,3,2,3) = matrix(5,5)/coe2;
      tensor(2,3,3,2) = matrix(5,5)/coe2;
      tensor(2,3,1,3) = matrix(5,6)/coe2;
      tensor(2,3,3,1) = matrix(5,6)/coe2;

      tensor(1,3,1,1) = matrix(6,1)/coe1;
      tensor(1,3,2,2) = matrix(6,2)/coe1;
      tensor(1,3,3,3) = matrix(6,3)/coe1;
      tensor(1,3,1,2) = matrix(6,4)/coe2;
      tensor(1,3,2,1) = matrix(6,4)/coe2;
      tensor(1,3,2,3) = matrix(6,5)/coe2;
      tensor(1,3,3,2) = matrix(6,5)/coe2;
      tensor(1,3,1,3) = matrix(6,6)/coe2;
      tensor(1,3,3,1) = matrix(6,6)/coe2;
      
      tensor(2,1,1,1) = matrix(4,1)/coe1;
      tensor(2,1,2,2) = matrix(4,2)/coe1;
      tensor(2,1,3,3) = matrix(4,3)/coe1;
      tensor(2,1,1,2) = matrix(4,4)/coe2;
      tensor(2,1,2,1) = matrix(4,4)/coe2;
      tensor(2,1,2,3) = matrix(4,5)/coe2;
      tensor(2,1,3,2) = matrix(4,5)/coe2;
      tensor(2,1,1,3) = matrix(4,6)/coe2;
      tensor(2,1,3,1) = matrix(4,6)/coe2;

      tensor(3,2,1,1) = matrix(5,1)/coe1;
      tensor(3,2,2,2) = matrix(5,2)/coe1;
      tensor(3,2,3,3) = matrix(5,3)/coe1;
      tensor(3,2,1,2) = matrix(5,4)/coe2;
      tensor(3,2,2,1) = matrix(5,4)/coe2;
      tensor(3,2,2,3) = matrix(5,5)/coe2;
      tensor(3,2,3,2) = matrix(5,5)/coe2;
      tensor(3,2,1,3) = matrix(5,6)/coe2;
      tensor(3,2,3,1) = matrix(5,6)/coe2;

      tensor(3,1,1,1) = matrix(6,1)/coe1;
      tensor(3,1,2,2) = matrix(6,2)/coe1;
      tensor(3,1,3,3) = matrix(6,3)/coe1;
      tensor(3,1,1,2) = matrix(6,4)/coe2;
      tensor(3,1,2,1) = matrix(6,4)/coe2;
      tensor(3,1,2,3) = matrix(6,5)/coe2;
      tensor(3,1,3,2) = matrix(6,5)/coe2;
      tensor(3,1,1,3) = matrix(6,6)/coe2;
      tensor(3,1,3,1) = matrix(6,6)/coe2;

end

function [ matrix ] = mat4_mat2( tensor,coe )

      matrix = zeros (6,6);
      if coe==1
          coe1=1;
          coe2=1;
      elseif coe==2
          coe1=2;
          coe2=4;
      end
      matrix(1,1)=tensor(1,1,1,1);
      matrix(1,2)=tensor(1,1,2,2);
      matrix(1,3)=tensor(1,1,3,3);
      matrix(1,4)=tensor(1,1,1,2)*coe1;
      matrix(1,5)=tensor(1,1,2,3)*coe1;
      matrix(1,6)=tensor(1,1,1,3)*coe1;

      matrix(2,1)=tensor(2,2,1,1);
      matrix(2,2)=tensor(2,2,2,2);
      matrix(2,3)=tensor(2,2,3,3);
      matrix(2,4)=tensor(2,2,1,2)*coe1;
      matrix(2,5)=tensor(2,2,2,3)*coe1;
      matrix(2,6)=tensor(2,2,1,3)*coe1;

      matrix(3,1)=tensor(3,3,1,1);
      matrix(3,2)=tensor(3,3,2,2);
      matrix(3,3)=tensor(3,3,3,3);
      matrix(3,4)=tensor(3,3,1,2)*coe1;
      matrix(3,5)=tensor(3,3,2,3)*coe1;
      matrix(3,6)=tensor(3,3,1,3)*coe1;

      matrix(4,1)=tensor(1,2,1,1)*coe1;
      matrix(4,2)=tensor(1,2,2,2)*coe1;
      matrix(4,3)=tensor(1,2,3,3)*coe1;
      matrix(4,4)=tensor(1,2,1,2)*coe2;
      matrix(4,5)=tensor(1,2,2,3)*coe2;
      matrix(4,6)=tensor(1,2,1,3)*coe2;

      matrix(5,1)=tensor(2,3,1,1)*coe1;
      matrix(5,2)=tensor(2,3,2,2)*coe1;
      matrix(5,3)=tensor(2,3,3,3)*coe1;
      matrix(5,4)=tensor(2,3,1,2)*coe2;
      matrix(5,5)=tensor(2,3,2,3)*coe2;
      matrix(5,6)=tensor(2,3,1,3)*coe2;

      matrix(6,1)=tensor(1,3,1,1)*coe1;
      matrix(6,2)=tensor(1,3,2,2)*coe1;
      matrix(6,3)=tensor(1,3,3,3)*coe1;
      matrix(6,4)=tensor(1,3,1,2)*coe2;
      matrix(6,5)=tensor(1,3,2,3)*coe2;
      matrix(6,6)=tensor(1,3,1,3)*coe2;
end

function [C] = Aijpq_Bpqkl(A,B)
C = zeros(3,3,3,3);
  for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for p=1:3
                    for q=1:3
                        C(i,j,k,l)=C(i,j,k,l)+A(i,j,p,q)*B(p,q,k,l);
                    end
                end
            end
        end
    end
  end

end

function [ C ] = Aijklpq_Bpq(A,B)
C = zeros(3,3,3,3);
  for i=1:3
      for j=1:3
        for k=1:3
            for l=1:3
                for p=1:3
                    for q=1:3
                        C(i,j,k,l)=C(i,j,k,l)+A(i,j,k,l,p,q)*B(p,q);
                    end
                end
            end
        end
    end
  end

end

function [matDz, matS] = matDO1(Omega)
global E0 nu0 a1 a2 a3 a4
b1 = (1+nu0)/E0/2;
b2 = nu0/E0;
E = eye(3);
E6 = eye(6);
trOmega = trace(Omega);
matS(1:3,1:3,1:3,1:3) = 0;
matS_2(1:6,1:6) = 0;
matDz(1:3,1:3,1:3,1:3) = 0;
matDz_2(1:6,1:6) = 0;
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                matS(i,j,k,l) = b1*(E(i,k)*E(j,l)+E(i,l)*E(j,k))-...
                    b2*E(i,j)*E(k,l)+2*a1*trOmega*E(i,j)*E(k,l)+...
                    0.5*a2*(E(i,k)*Omega(j,l)+E(i,l)*Omega(j,k)+...
                    Omega(i,k)*E(j,l)+Omega(i,l)*E(j,k))+...
                    a3*(E(i,j)*Omega(k,l)+Omega(i,j)*E(k,l))+...
                    a4*trOmega*(E(i,k)*E(j,l)+E(i,l)*E(j,k));
            end
        end
    end
end    
matS_2 = mat4_mat2(matS,2);
matDz_2 = matS_2\E6;  % Give inverse matrix
matDz = mat2_mat4(matDz_2,1);

end

function [ P_1 ] = matP_1(sigmaT)
%   P_1 projection tensor
%   P_1=H(sigmaT)-H(-sigmaT)

P_1(1:3,1:3,1:3,1:3)=0;

n1=[1;0;0];
n2=[0;1;0];
n3=[0;0;1];
[V,D] = eig(sigmaT);
n1=V(:,1);
n2=V(:,2);
n3=V(:,3);
P_1=(heaviside(D(1,1))-heaviside(-D(1,1)))*ni_nj_nk_nl(n1)+...
    (heaviside(D(2,2))-heaviside(-D(2,2)))*ni_nj_nk_nl(n2)+...
    (heaviside(D(3,3))-heaviside(-D(3,3)))*ni_nj_nk_nl(n3);

end

function [ P_2 ] = matP_2(sigmaT)

P_2(1:3,1:3,1:3,1:3)=0;

n1=[1;0;0];
n2=[0;1;0];
n3=[0;0;1];
[V,D] = eig(sigmaT);
n1=V(:,1);
n2=V(:,2);
n3=V(:,3);
s(1:3)=0;
s(1)=D(1,1)-max(D(1,1),max(D(2,2),D(3,3)));
s(2)=D(2,2)-max(D(1,1),max(D(2,2),D(3,3)));
s(3)=D(3,3)-max(D(1,1),max(D(2,2),D(3,3)));
for ii=1:3
    if s(ii)==0
        s(ii)=0;
    else
        s(ii)=heaviside(-s(ii));
    end
end

P_2=s(1)*ni_nj_nk_nl(n1)+...
    s(2)*ni_nj_nk_nl(n2)+...
    s(3)*ni_nj_nk_nl(n3);

end

function [ M ] = ni_nj_nk_nl( n )

M(1:3,1:3,1:3,1:3)=0;

for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                M(i,j,k,l) = n(i)*n(j)*n(k)*n(l);
            end
        end
    end
end 
end

function [ fdm,sigm,Omegam,depsidm ] = fd_lam( lambdam, Omega, sigmaT, epsT, epsid, Pars,deps)
% Iteration solving for lambdan

a1=Pars(1);
a2=Pars(2);
a3=Pars(3);
a4=Pars(4);
C0=Pars(5);
C1=Pars(6);
alpha=Pars(7);

E=eye(3);
Yd1 = a1*(trace(sigmaT))^2*E+a2*sigmaT*sigmaT+a3*trace(sigmaT)*sigmaT+a4*trace(sigmaT*sigmaT)*E;
dY_dsig = dY_dsigf(sigmaT);
P_1 = matP_1(sigmaT);
P_2 = matP_2(sigmaT);
F1ij = Aijkl_Bkl(P_1,Yd1)-1/3*Aij_Bij(Aijkl_Bkl(P_1,Yd1),E)*E;
F2ij = Aijkl_Bkl(P_2,Yd1);
df_dOmega = -C1*E;
F2F2=Aij_Bij(F2ij,F2ij);
if F2F2==0
    dg_dY=zeros(3,3);
else    
    dg_dY = Aijkl_Bij(P_2,F2ij)/sqrt(2*Aij_Bij(F2ij,F2ij));
end
dY_dsig = dY_dsigf(sigmaT);
P_3 = P_1-1/3*Aij_Bkl(E,Aijkl_Bij(P_1,E));
df_dY = Aijkl_Bij(P_3,F1ij)/sqrt(2*Aij_Bij(F1ij,F1ij))-alpha*Aijkl_Bij(P_1,E);
% updating depsid
depsidm = lambdam*Aijkl_Bij(dY_dsig,df_dY);
% depsidn = lambdan*Aijkl_Bij(dY_dsig,df_dY);
% updating Omega
dOmegam = lambdam*dg_dY;
% dOmegan = lambdan*dg_dY;
Omegam = Omega+dOmegam;
% Omegan = Omega+dOmegan;
[matD,S] = matDO1(Omega);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for p=1:3
                    for q=1:3
                        dS_dOmega(i,j,k,l,p,q)=2*a1*E(i,j)*E(k,l)*E(p,q)+...
                            1/4*a2*(E(i,k)*(E(p,j)*E(q,l)+E(p,l)*E(q,j))+...
                        E(i,l)*(E(p,j)*E(q,k)+E(p,k)*E(q,j))+...
                        E(j,l)*(E(i,p)*E(q,k)+E(i,q)*E(p,k))+...
                        E(j,k)*(E(i,p)*E(q,l)+E(i,q)*E(p,l)))+...
                        1/2*a3*(E(i,j)*(E(k,p)*E(l,q)+E(k,q)*E(l,p))+E(k,l)*(E(i,p)*E(j,q)+E(i,q)*E(j,p)))+...
                        a4*(E(i,k)*E(j,l)+E(i,l)*E(j,k))*E(p,q);
                    end
                end
            end
        end
    end
end
temp1 = -invmat4(Aijpq_Bpqkl(S,S));
dmatD_dOmega=zeros(3,3,3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for m=1:3
                    for n=1:3
                        for p=1:3
                            for q=1:3
    dmatD_dOmega(i,j,k,l,m,n)=dmatD_dOmega(i,j,k,l,m,n)+temp1(i,j,p,q)*dS_dOmega(p,q,k,l,m,n);
                            end
                        end
                    end
                end
            end
        end
    end
end

matDem = matD + lambdam*Aijklpq_Bpq(dmatD_dOmega,dg_dY);
epsm = epsT+deps-epsid-depsidm;
sigm = Aijkl_Bkl(matDem,epsm);

[fdm,Yd1] = fdDP(sigm,zeros(3,3),Omegam,Pars);
end

function [ fdm,epsm,Omegam,depsidm ] = fd_lam2( lambdam, Omega, sigmaT, epsid, Pars,dsig)
%    Iteration solving for lambdan
%

a1=Pars(1);
a2=Pars(2);
a3=Pars(3);
a4=Pars(4);
C0=Pars(5);
C1=Pars(6);
alpha=Pars(7);

E=eye(3);
Yd1 = a1*(trace(sigmaT))^2*E+a2*sigmaT*sigmaT+a3*trace(sigmaT)*sigmaT+a4*trace(sigmaT*sigmaT)*E;
dY_dsig = dY_dsigf(sigmaT);
P_1 = matP_1(sigmaT);
P_2 = matP_2(sigmaT);
F1ij = Aijkl_Bkl(P_1,Yd1)-1/3*Aij_Bij(Aijkl_Bkl(P_1,Yd1),E)*E;
F2ij = Aijkl_Bkl(P_2,Yd1);
df_dOmega = -C1*E;
F2F2=Aij_Bij(F2ij,F2ij);
if F2F2==0
    dg_dY=zeros(3,3);
else    
    dg_dY = Aijkl_Bij(P_2,F2ij)/sqrt(2*Aij_Bij(F2ij,F2ij));
end
dY_dsig = dY_dsigf(sigmaT);
P_3 = P_1-1/3*Aij_Bkl(E,Aijkl_Bij(P_1,E));
df_dY = Aijkl_Bij(P_3,F1ij)/sqrt(2*Aij_Bij(F1ij,F1ij))-alpha*Aijkl_Bij(P_1,E);
% updating depsid
depsidm = lambdam*Aijkl_Bij(dY_dsig,df_dY);
% depsidn = lambdan*Aijkl_Bij(dY_dsig,df_dY);
% updating Omega
dOmegam = lambdam*dg_dY;
% dOmegan = lambdan*dg_dY;
Omegam = Omega+dOmegam;
% Omegan = Omega+dOmegan;
[matD,S] = matDO1(Omega);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for p=1:3
                    for q=1:3
                        dS_dOmega(i,j,k,l,p,q)=2*a1*E(i,j)*E(k,l)*E(p,q)+...
                            1/4*a2*(E(i,k)*(E(p,j)*E(q,l)+E(p,l)*E(q,j))+...
                        E(i,l)*(E(p,j)*E(q,k)+E(p,k)*E(q,j))+...
                        E(j,l)*(E(i,p)*E(q,k)+E(i,q)*E(p,k))+...
                        E(j,k)*(E(i,p)*E(q,l)+E(i,q)*E(p,l)))+...
                        1/2*a3*(E(i,j)*(E(k,p)*E(l,q)+E(k,q)*E(l,p))+E(k,l)*(E(i,p)*E(j,q)+E(i,q)*E(j,p)))+...
                        a4*(E(i,k)*E(j,l)+E(i,l)*E(j,k))*E(p,q);
                    end
                end
            end
        end
    end
end

Sm=zeros(3,3,3,3);
Sm=S+Aijklpq_Bpq(dS_dOmega,dOmegam);

sigm = sigmaT+dsig;
epsm = Aijkl_Bkl(Sm,sigm)+epsid+depsidm;
[fdm,Yd1] = fdDP(sigm,zeros(3,3),Omegam,Pars);

end
