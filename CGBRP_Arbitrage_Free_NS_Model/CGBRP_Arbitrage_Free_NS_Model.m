%data input
xdata=xlsread('Yields.xlsx');
colnumber = 2:1:12;
y=xdata(1:end-3,colnumber)/1200;
t=size(y,1);
k=size(y,2);


%Define yield structure used in the estimation

yname =  ['m01';'m02';'m03';'m06';'m09';'y01';'y02';'y03';'y05';'y07';'y10']; 

maturity = [ 1;    2;    3;    6;    9;   12;   24;   36;   60;   84;  120];
yn=size(yname,1);


%mcmc
irep1=0;
irep2=100;
irep=irep1+irep2;

%store
mu_l_q_=[];
omega_=zeros(irep2,3,3);;
phi=zeros(irep2,3,3);
lamda_=[];
mu_=[];
x_=zeros(irep2,t,3);;
h_=[];
cm=zeros(11,1);
dosmo=1;

%initial value
mu_l_q=0;
mu_q=[mu_l_q,0,0]';
omega= 10^(-5)* [0.0078 -0.0023 -0.0021; -0.0023 0.0277 0.0010; -0.0021 0.001 0.1069];
phi=[0.9863 0.0231 -0.0055;-0.0122 0.925 0.0219;0.1784 0.2057 0.7450];
lamda=0.077;
mu=10^(-3)*[0.1281 -0.0495 -0.7853]';
sigma=1.5951*10^(-4);
h=1/sigma^(2);
for i=1:irep
i
%xt part
[A BB]=coeff(yn,maturity,lamda,mu_q,omega);
theta=inv(eye(3)-phi)*mu;
%for demean
for j=1:yn
demean(j,1)=A(j,1)+BB(j,:)*theta;
end
for j=1:yn
ydemean(:,j)=y(:,j)-ones(size(y,1),1)*demean(j,1);
end
sigma2=1/h;
[mx xt Ptt] = ssp2(ydemean,BB,cm,phi,omega,sigma2,dosmo);
Theta=[];
for j=1:t
Theta=[Theta theta];
end
x=xt+Theta;


%mu,phi,omega part
X=[ones(size(y,1)-1,1) x(:,1:end-1)'];
Z=[x(:,2:end)]';
phi_hat=inv(X'*X)*(X'*Z);
H_hat=(Z-X*phi_hat)'*(Z-X*phi_hat);
tau=size(y,1)-4;
K=3;


omega = inv(wish(inv(H_hat),tau));% Draw SIGMA


v_phi=kron(H_hat,inv(X'*X))/tau;
mu_phi=reshape(phi_hat,12,1);
Phi=mvnrnd(mu_phi,v_phi); 

%reshape
Phi_matrix=reshape(Phi,4,3);
mu=Phi_matrix(1,:)';
phi=Phi_matrix([2 3 4],:)';

%lamda,mu_l_q,sigma2 part

%h part
[A BB]=coeff(yn,maturity,lamda,mu_q,omega);
s2=0;
for j=1:size(y,1)
s2=s2+(y(j,:)'- BB*x(:,j) -A)'*(y(j,:)'- BB*x(:,j) -A);
end

v=size(y,1)*size(y,2);
%h=gamrnd(v/s2,v);
h=gamm_rnd(1,1,0.5*v,0.5*s2);

%lamda part
lamda_draw=lamda+0.0001*randn(1,1);
[A_lamda_draw BB_lamda_draw]=coeff(yn,maturity,lamda_draw,mu_q,omega);
[l] = likelihood(y,x,A,BB,h);
[l_lamda_draw] = likelihood(y,x,A_lamda_draw,BB_lamda_draw,h);
laccprob = l_lamda_draw-l;
   if log(rand)<laccprob
       lamda=lamda_draw;
   end 

%mu_l_q part
mu_l_q_draw=mu_l_q+0.0001*randn(1,1);
mu_q_draw=[mu_l_q_draw,0,0]';

[A BB]=coeff(yn,maturity,lamda,mu_q,omega);
[A_mu_l_q_draw BB_mu_l_q_draw]=coeff(yn,maturity,lamda,mu_q_draw,omega);

[l] = likelihood(y,x,A,BB,h);
[l_mu_l_q_draw] = likelihood(y,x,A_mu_l_q_draw,BB_mu_l_q_draw,h);

laccprob = l_mu_l_q_draw-l; 

%accept candidate draw with log prob = laccprob, else keep old draw
   if log(rand)<laccprob
       mu_l_q=mu_l_q_draw;
   end 
mu_q=[mu_l_q,0,0]';


    if i>irep1
        %after discarding burnin, store all draws
        phi_(i-irep1,:,:)=phi;
        x_(i-irep1,:,:)=x';
       omega_(i-irep1,:,:)=omega;
    end
end


%Risk premium calculation

for i=1:t
for j=1:3
x50(i,j)=median(x_(i,j,:));
end
end
x50=x50';

 for i=1:3
  mu_post(i,:)=mean(mu_(i,:));
  end

 for i=1:3
 for j=1:3
 phi_post(i,j)=mean(phi_(i,j,:));
 end
 end

 for i=1:3
 for j=1:3
 omega_post(i,j)=mean(omega_(i,j,:));
 end
 end

yield_risk_neutral=(BB*x50)';

for i=1:t
yield_risk_neutral(i,:)=yield_risk_neutral(i,:)+A';
end

[A_ex BB_ex] = coeff_risk(yn,maturity,lamda,mu_post,phi_post,omega_post);
yield_ex=(BB_ex*x50)';

for i=1:t
yield_ex(i,:)=yield_ex(i,:)+A_ex';
end

risk_premium_10=yield_risk_neutral(:,end)-yield_ex(:,end);