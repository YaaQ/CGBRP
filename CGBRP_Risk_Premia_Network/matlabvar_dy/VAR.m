function [b, vc] = VAR(data, lag)

    p=lag ;

    y=data ;
    ylag = lagmatrix(y,1:p) ;
    ylag(any(isnan(ylag),2),:) = [];
    
    n=size(data,2) ;

    xmat = [ones(length(ylag),1) ylag];
    ymat = trimr(y,p,0);

    bar = xmat\ymat;
    nconst = 2;
    mu = bar(1:nconst,:);

    v   = ymat - xmat*bar;
	vc = v'*v/(size(v,1)-(n+1)*lag);
    
    b = transpose(bar) ;



end