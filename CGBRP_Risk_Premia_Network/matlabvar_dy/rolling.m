
data_raw = xlsread('C:\Users\demirermert\Desktop\connectedness\Matlab\data.xlsx') ;

data=log(data_raw) ;

obs = size(data,1) ;

n_country = size(data,2) ;

roll = 100 ;

lag = 3 ;

h = 10 ;

net_=zeros(n_country, obs-roll+1-lag) ; 
from_=zeros(n_country, obs-roll+1-lag) ;
to_=zeros(n_country, obs-roll+1-lag) ;
index_=zeros(obs-roll-lag+1,1) ;
table_=zeros(n_country, n_country, obs-roll-lag+1) ;

for i=1:(obs-roll-lag+1)
    
    data_roll=data(i:(i+roll+lag-1),:) ;
    
    [b, vc] = VAR(data ,lag) ;
    
    [to, from, net, total, con] = index(b, vc, h, lag) ;
    
    to_(:,i) = to ;
    from_(:,i) = from ;
    net_(:,i) = net ;
    index_(i,1) = total ;
    table_(:,:,i) = con ;
    
end

