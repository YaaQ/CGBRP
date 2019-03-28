
addpath(genpath('/Users/local/Desktop/matlab/matlabvar_dy'))
load('risk.mat')
risk=[ca_risk(:,end) ge_risk(:,end) uk_risk(:,end) usa_risk(:,end)];
to_=[];
from_=[];
net_=[];
total_=[];
con_=[];

lag=1;
h=12;
for i=1:202-120
data=risk(i:120+i,:);
[b, vc] = VAR(data, lag)
[to, from, net, total, con] = index(b, vc, h ,lag)
to_=[to_ to];
from_=[from_ from];
net_=[net_ net];
total_=[total_ total];
con_=[con_];
end



%-----------------------
[a,b,c,d,e]=textread('risk_premium_time.txt','%d %1s %d %1s %d');
clear b d
ttc=datenum([a c e]);
dateaxis('x',2)

%------------------
plot(ttc(121:end),net_(1,:))
dateaxis('x',2)


subplot(4,1,1)
plot(ttc,risk(:,2),'-k')
 dateaxis('x',2)
 subplot(4,1,2)
plot(ttc,risk(:,4),'-r')
dateaxis('x',2)
 subplot(4,1,3)
plot(ttc,risk(:,3),'-y')
dateaxis('x',2)
subplot(4,1,4)
plot(ttc,risk(:,1),'-b')
dateaxis('x',2)




plot(maturity_germany,germany(end,:),'-k');
Hold on
plot(maturity_usa,usa(end,:),'-r');
plot(maturity_uk,uk(end,:),'-y');
plot(maturity_c,canada(end,:),'-b');