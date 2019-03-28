function [to, from, net, total, con] = index(b, vc, h ,lag)

    n=size(vc,2) ;

    if size(b,2)==lag*n
    
        b = [zeros(size(b,1),1), b] ;
        
    end
    
    for k=0:(n-1)
        
        if k>0
            
            vc_temp = vc ;
            b_temp = b ;
            b_temp = b(1,:) ;
            b = [b(2:n,:); b_temp] ;
             
       

            for i=0:(lag-1)

                b_shift = b(:,(n*i+2):(n*(i+1)+1)) ;
                b_temp = b_shift(:,1) ;
                b_shift = [b_shift(:,2:n) b_temp] ;
                b(:,(n*i+2):(n*(i+1)+1)) = b_shift ;

            end

            for i=0:(n-1)

                for j=0:(n-1)

                    u=i ;
                    g=j ;

                    if u==0                  
                        u=n ;                  
                    end

                    if g==0                  
                        g=n ;                  
                    end  

                
                    vc_temp(u,g) = vc((i+1),(j+1)) ;

                end

            end

            vc = vc_temp ;

        end

        si = zeros(n, (h+lag-1)*n) ; 

        m = eye(n) ;

        si(:,((lag-1)*n+1):(lag*n)) = m ;

        for i=lag:(h +lag-2)

            for p=0:(lag-1) 

                si(:,(i*n+1):((i+1)*n)) = b(:,(n*p+2):(n*(p+1)+1)) * si(:,((i-p-1)*n+1):((i-p)*n)) + si(:,(i*n+1):((i+1)*n))  ;         

            end

        end

        si = si(:,(((lag-1)*n)+1):size(si,2)) ;
        

        s = chol(vc,'lower') ;

        irf = zeros(n, n*h) ;

        for i=0:(h-1)

           irf(:,(i*n+1):((i+1)*n)) = si(:,(i*n+1):((i+1)*n)) * s ;      

        end


        irf2 = irf .* irf ;
        

        if k<1
      
          sum_ = sum(irf2(:,(1:n)),2) * h ;

        end   

        dec=zeros(n, n) ;

        for i=0:(h-1)

          dec=dec+irf2(:,(i*n+1):((i+1)*n)) ;

        end

        if k>0 

          temp = dec((n-k+1):n,:) ;

          dec=[temp ; dec(1:(n-k),:)]   ;

        end

        con(:,(k+1)) = dec(:,1) ;
    
    end
    


    for l=0:(n-1)

    con((l+1),:) = con((l+1),:)/sum_(l+1,1)  ;
    con((l+1),:) = con((l+1),:)./sum(con((l+1),:),2) ;

    end
	
	con_dif = con - transpose(con) ;
    
    r = diag(con) ;
    
    t = sum(con,1) ;

	from=((sum(con,2))-diag(con))*100  ;
	to=(transpose(sum(con,1))-diag(con))*100	 ;
	net=to-from ;
	total=sum(to,1)/n ;


end

