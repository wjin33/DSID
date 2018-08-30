function Kachanov()

clear all; close all;
x0=0;y0=0;
dx=50;dy=50;
density(1,1:4)=0.45:0.05:0.6;
YonModu(1,1:4)=0;
SheModu(1,1:4)=0;
l=1;
for kkk=4
    jointNo=round((dx-x0)*(dy-y0)*density(1,kkk));
    CrackInfo=zeros(4,jointNo);
    i=1;
while i<=jointNo
    xc=0.4+0.2*round(rand*246);
    yc=50*rand;
    y1=yc-l;y2=yc+l;
    count=0;
    if y1 >= (y0+0.4)  && y2 <= (y0+dy-0.4)
        if i>1
            count_full=0;
            count=0;
            for j=1:i-1
                if xc==CrackInfo(1,j)
                    count_full=count_full+1;
                    if yc+1.5 <= CrackInfo(2,j)-1 || yc-1.5 >= CrackInfo(2,j)+1
                        count=count+1;
                    else
                        break
                    end
                end
            end
            if count_full==count
                CrackInfo(:,i)=[xc,yc,pi/2,l];
                i=i+1;
            end     
        else
           CrackInfo(:,i)=[xc,yc,pi/2,l];
           i=i+1;
        end
    end
end
% end
% for kkk=1:4
%     l=(kkk+8)/12;
%     CrackInfo(4,:)=l;
%Calculating the Lambda matrix
Lammda=zeros(2*jointNo,2*jointNo);
 for i=1:jointNo
     llocal=-CrackInfo(4,i):0.005:CrackInfo(4,i);
     [~,PointNo]=size(llocal);
     for j=1:jointNo      
         if i~=j
             theta=CrackInfo(3,i)-CrackInfo(3,j);
             alpha=CrackInfo(3,j);
             xlocalij=(CrackInfo(1,i)-CrackInfo(1,j))*cos(alpha)+...   %%%%%modify 0
             (CrackInfo(2,i)-CrackInfo(2,j))*sin(alpha)+llocal.*cos(theta); %The influence of crack j on i;
             ylocalij=-(CrackInfo(1,i)-CrackInfo(1,j))*sin(alpha)+...
             (CrackInfo(2,i)-CrackInfo(2,j))*cos(alpha)+llocal.*sin(theta);   %%%%%modify 0
             [sigma_xx,sigma_xy,sigma_yy,tau_xx,tau_xy,tau_yy]=addtional_traction(CrackInfo(4,j),xlocalij,ylocalij);
             sigma_yy=real(sigma_yy);
             sigma_xx=real(sigma_xx);
             sigma_xy=real(sigma_xy);
             tau_yy=real(tau_yy);
             tau_xx=real(tau_xx);
             tau_xy=real(tau_xy);
             if theta>=pi/2
                ni_local=[sin(theta);-cos(theta)];
                gammai_local=[cos(theta);sin(theta)];
             else
                ni_local=[-sin(theta);cos(theta)];
                gammai_local=[cos(theta);sin(theta)];
             end
                temp=zeros(4,PointNo);
             for k=1:PointNo
                 stress_sigma=[sigma_xx(1,k),sigma_xy(1,k);sigma_xy(1,k),sigma_yy(1,k)];
                 stress_tau=[tau_xx(1,k),tau_xy(1,k);tau_xy(1,k),tau_yy(1,k)];
                 temp(1,k)=ni_local'*stress_sigma*ni_local;
                 temp(2,k)=ni_local'*stress_tau*ni_local;
                 temp(3,k)=gammai_local'*stress_sigma*ni_local;
                 temp(4,k)=gammai_local'*stress_tau*ni_local;
             end
             Lammda(i*2-1,j*2-1) = -mean(temp(1,:));
             Lammda(i*2-1,j*2)   = -mean(temp(2,:));
             Lammda(i*2,j*2-1)   = -mean(temp(3,:));
             Lammda(i*2,j*2)     = -mean(temp(4,:));
         else
             Lammda(i*2-1,j*2-1) = 1;
             Lammda(i*2-1,j*2)   = 0;
             Lammda(i*2,j*2-1)   = 0;
             Lammda(i*2,j*2)     = 1;
         end
     end
 end
 %%%%%%
 clear temp llocal xlocalij ylocalij sigma_xx sigma_xy sigma_yy tau_xx tau_xy tau_yy stress_sigma stress_tau
 %%%%%%%%%%%
 sigma=[1,0;0,0];
 p_tau_infinite=zeros(2*jointNo,1);
 for i=1:jointNo
     gamma=[cos(CrackInfo(3,i)); sin(CrackInfo(3,i))];
     n=[-sin(CrackInfo(3,i));cos(CrackInfo(3,i))];
     p_tau_infinite(i*2-1,1)=  n'*sigma*n;
     p_tau_infinite(i*2,1)= gamma'*sigma*n; 
 end
p_tau_local = Lammda\p_tau_infinite;
% clear Lammda p_tau_infinite gamma n

E0 = 4.6E10;  % unit=Pa
nu0 = 0.186;

crack_strain=zeros(2,2);
for i=1:jointNo
    if CrackInfo(3,i)>=pi/2
       ni=[sin(CrackInfo(3,i));-cos(CrackInfo(3,i))];
       else
       ni=[-sin(CrackInfo(3,i));cos(CrackInfo(3,i))];
    end
    bi=[p_tau_local(2*i-1,1); p_tau_local(2*i,1)];
    Transform=[sin(CrackInfo(3,i))  cos(CrackInfo(3,i))
               cos(CrackInfo(3,i)) -sin(CrackInfo(3,i));];
    bi=(Transform*bi).*(pi*CrackInfo(4,i)/E0);  
    crack_strain(1,1)=crack_strain(1,1)+(ni(1)*bi(1)+bi(1)*ni(1))*CrackInfo(4,i)*2/(2*(dy-y0)*(dx-x0));
    crack_strain(1,2)=crack_strain(1,2)+(ni(1)*bi(2)+bi(1)*ni(2))*CrackInfo(4,i)*2/(2*(dy-y0)*(dx-x0));
    crack_strain(2,1)=crack_strain(2,1)+(ni(2)*bi(1)+bi(2)*ni(1))*CrackInfo(4,i)*2/(2*(dy-y0)*(dx-x0));
    crack_strain(2,2)=crack_strain(2,2)+(ni(2)*bi(2)+bi(2)*ni(2))*CrackInfo(4,i)*2/(2*(dy-y0)*(dx-x0));
end
elastic_strain=zeros(2,2);
elastic_strain(1,1)=(sigma(1,1)-nu0*sigma(2,2))/E0;
elastic_strain(2,2)=(sigma(2,2)-nu0*sigma(1,1))/E0;
elastic_strain(1,2)=sigma(1,2)*(1+nu0)/E0;

total_strain_vector=[elastic_strain(1,1)+crack_strain(1,1);
                     elastic_strain(2,2)+crack_strain(2,2);
                     elastic_strain(1,2)+crack_strain(1,2)];
stress_vector=[sigma(1,1) sigma(2,2) sigma(1,2)];
Ee1=stress_vector(1,1)/total_strain_vector(1,1);

%%%%%%%%%%%
 sigma=[0,1;1,0];
 p_tau_infinite=zeros(2*jointNo,1);
 for i=1:jointNo
     gamma=[cos(CrackInfo(3,i)); sin(CrackInfo(3,i))];
     n=[-sin(CrackInfo(3,i));cos(CrackInfo(3,i))];
     p_tau_infinite(i*2-1,1)=  n'*sigma*n;
     p_tau_infinite(i*2,1)= gamma'*sigma*n; 
 end
p_tau_local = Lammda\p_tau_infinite;
clear Lammda p_tau_infinite gamma n

E0 = 4.6E10;  % unit=Pa
nu0 = 0.186;

crack_strain=zeros(2,2);
for i=1:jointNo
    if CrackInfo(3,i)>=pi/2
       ni=[sin(CrackInfo(3,i));-cos(CrackInfo(3,i))];
       else
       ni=[-sin(CrackInfo(3,i));cos(CrackInfo(3,i))];
    end
    bi=[p_tau_local(2*i-1,1); p_tau_local(2*i,1)];
    Transform=[sin(CrackInfo(3,i))  cos(CrackInfo(3,i))
               cos(CrackInfo(3,i)) -sin(CrackInfo(3,i));];
    bi=(Transform*bi).*(pi*CrackInfo(4,i)/E0);  
    crack_strain(1,1)=crack_strain(1,1)+(ni(1)*bi(1)+bi(1)*ni(1))*CrackInfo(4,i)*2/(2*(dy-y0)*(dx-x0));
    crack_strain(1,2)=crack_strain(1,2)+(ni(1)*bi(2)+bi(1)*ni(2))*CrackInfo(4,i)*2/(2*(dy-y0)*(dx-x0));
    crack_strain(2,1)=crack_strain(2,1)+(ni(2)*bi(1)+bi(2)*ni(1))*CrackInfo(4,i)*2/(2*(dy-y0)*(dx-x0));
    crack_strain(2,2)=crack_strain(2,2)+(ni(2)*bi(2)+bi(2)*ni(2))*CrackInfo(4,i)*2/(2*(dy-y0)*(dx-x0));
end
elastic_strain=zeros(2,2);
elastic_strain(1,1)=(sigma(1,1)-nu0*sigma(2,2))/E0;
elastic_strain(2,2)=(sigma(2,2)-nu0*sigma(1,1))/E0;
elastic_strain(1,2)=sigma(1,2)*(1+nu0)/E0;

total_strain_vector=[elastic_strain(1,1)+crack_strain(1,1);
                     elastic_strain(2,2)+crack_strain(2,2);
                     elastic_strain(1,2)+crack_strain(1,2)];
stress_vector=[sigma(1,1) sigma(2,2) sigma(1,2)];
Ge12=stress_vector(1,3)/(2*total_strain_vector(3,1));


    YonModu(1,kkk)=Ee1/E0;
    G0=E0/(2*(1+nu0));
    SheModu(1,kkk)=Ge12/G0;
%   kkk=kkk+1;
    save case_random_init_addtion4 density YonModu SheModu
end
save case_random_init_addtion4 density YonModu SheModu

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sigma_xx,sigma_xy,sigma_yy,tau_xx,tau_xy,tau_yy]=addtional_traction(l,x,y)
% clear all; close all;
% l=0.5;
% x=-1.6:0.01:-0.8;
% y=zeros(size(x));

alpha=(x-l).^2+y.^2;
beta=2*(x.^2+y.^2-l^2);
gamma=(x+l).^2+y.^2;
delta=beta+2*sqrt(alpha.*gamma);


I1=4*l^3*(sqrt(gamma)-sqrt(alpha))./(sqrt(delta).*(sqrt(alpha)+sqrt(gamma)+sqrt(delta)).^2);
I2=4*l^2./(sqrt(delta).*(sqrt(alpha)+sqrt(gamma)+sqrt(delta)));
I3=2*l^3*(sqrt(gamma)-sqrt(alpha))./(sqrt(alpha.*gamma).*(delta).^(3/2));
I4=2*l^2*(sqrt(alpha)+sqrt(gamma))./(sqrt(alpha.*gamma).*(delta).^(3/2));
I5=(l^3/2)*(3*sqrt(alpha.*gamma).*(sqrt(alpha)+sqrt(gamma)).^2.*(sqrt(gamma)-sqrt(alpha))...
    +delta.*((gamma).^(3/2)-(alpha).^(3/2)))./((alpha.*gamma).^(3/2).*(delta).^(5/2));
I6=(l^2/2).*(3*sqrt(alpha.*gamma).*(sqrt(alpha)+sqrt(gamma)).^3+...
    +delta.*((alpha).^(3/2)+(gamma).^(3/2)))./((alpha.*gamma).^(3/2).*(delta).^(5/2));



sigma_xx=I2-8*y.^2.*I4+8*y.^4.*I6;
sigma_xy=2*(-y.*I3+x.*y.*I4+4*y.^3.*I5-4.*x.*y.^3.*I6);
sigma_yy=I2+4*y.^2.*I4-8*y.^4.*I6;


tau_xx=2*(3*y.*I3-3.*x.*y.*I4-4*y.^3.*I5+4.*x.*y.^3.*I6);
tau_xy=I2-8*y.^2.*I4+8.*y.^4.*I6;
tau_yy=2*(-y.*I3+x.*y.*I4+4.*y.^3.*I5-4*x.*y.^3.*I6);


% plot(x,sigma_xx)
% figure
% plot(x,sigma_xy)
% figure
% plot(x,sigma_yy)
% figure
% plot(x,tau_xx)
% figure
% plot(x,tau_xy)
% figure
% plot(x,tau_yy)


end








