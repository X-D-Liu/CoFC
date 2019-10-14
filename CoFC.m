clear all;close all;clc

load ecoli.mat

points     = Y(:,[1:7]);
points_n   = size(points,1);
points_dim = size(points,2);

X=[]; 
[points_n,points_dim] = size(points); 

label=Y(:,8);
cluster_n = length(unique(Y(:,8)));
P = 2;
alpha = 0.02*ones(P,P);

X{1} = Y(:,[1,3,4,5,7]); 
X{2} = Y(:,[2,4,6]);
S{1} = size(X{1},2);
S{2} = size(X{2},2);

rng(10)
[clust_cen{1},u{1},obj1] = fcm(X{1},cluster_n); 
[clust_cen{2},u{2},obj2] = fcm(X{2},cluster_n);


for k = 1:cluster_n 
    D = zeros(1,cluster_n); 
    for l = 1:cluster_n 
        D(l) = sum((u{1}(l,:) - u{2}(l,:)).^2); 
    end 
    [minvalue,index] = min(D); 
    u{3}(index,:) = u{2}(l,:); 
    clust_cen{3}(index,:) = clust_cen{2}(l,:); 
end 
u{2} = u{3}; 
clust_cen{2} = clust_cen{3}; 

u_ref1 = u{1}; 
u_ref2 = u{2};


SSIGMA = 0; 
for itr = 1:50  
    itr
    
    
    for kk = 1:P
    
        for i = 1:points_n 
            for k = 1:cluster_n 
                temp1 = 0; 
                temp2 = 0;
                
                for jj = 1:P 
                    if jj ~= kk 
                        temp1 = temp1 + alpha(kk,jj)*u{jj}(k,i); 
                        temp2 = temp2 + alpha(kk,jj); 
                    end 
                end 
                
                temp3 = 0; 
                for l = 1:cluster_n 
                    temp3 = temp3 + (GetDistance(X{kk}(i,:),clust_cen{kk}(k,:)).^2)/(GetDistance(X{kk}(i,:),clust_cen{kk}(l,:)).^2); 
                end 
                
                temp4 = 0; 
                for l = 1:cluster_n 
                    temp = 0; 
                    for jj = 1:P 
                        if jj ~= kk 
                            temp = temp + alpha(kk,jj)*u{jj}(l,i); 
                        end 
                    end 
                     
                    temp4 = temp4 + temp/(1 + temp2); 
                end 
                
                u{kk}(k,i) = temp1/(1 + temp2) + (1 - temp4)/temp3; 
            end 
            
            u{kk}(:,i) = u{kk}(:,i)/sum(u{kk}(:,i));    
        end 
        
        
        for k = 1:cluster_n 
            for s = 1:S{kk} 
                temp1 = 0; 
                temp2 = 0; 
                temp3 = 0; 
                temp4 = 0; 
                for i = 1:points_n 
                    temp1 = temp1 + u{kk}(k,i)^2*X{kk}(i,s); 
                    temp2 = temp2 + u{kk}(k,i)^2; 
                    for jj = 1:P 
                        if jj ~= kk 
                            temp3 = temp3 + alpha(kk,jj)*(u{kk}(k,i) - u{jj}(k,i))^2*X{kk}(i,s); 
                            temp4 = temp4 + alpha(kk,jj)*(u{kk}(k,i) - u{jj}(k,i))^2; 
                        end 
                    end 
                end 
                clust_cen{kk}(k,s) = (temp1 + temp3)/(temp2 + temp4); 
            end 
        end 
            
    end
    
    SIGMA(itr) = sum(sum(abs(u{1} - u{2})))/(points_n*cluster_n); 
    if abs(SSIGMA - SIGMA(itr)) < 0.001 
        break; 
    else 
        SSIGMA = SIGMA(itr); 
    end 
    
    
end



[max_u,index] = max(u{1}) 
u=u{1}; 
u0=zeros(cluster_n,points_n); 
    for k=1:143 
        u0(1,k)=1; 
    end 
    for k=144:259 
        u0(2,k)=1; 
    end 
    for k=260:311 
        u0(3,k)=1; 
    end 
    for k=312:336 
        u0(4,k)=1; 
    end 

clust=[];
u_new=zeros(cluster_n,points_n); 
for k=1:size(u_new,1) 
    for l=1:size(u_new,2) 
        if u(k,l)==max(u(:,l)) 
           u_new(k,l)=1; 
           [num idx]=max(u(:,l));
           clust=[clust;idx];
        end 
    end 
end 



AR=1-ErrorRate(label,clust,cluster_n)/points_n
