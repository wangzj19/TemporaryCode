function result = PIO()
%%####A new bio-inspired swarm intelligence optimizer ?C pigeon inspired 
%%####optimization (PIO) is presented by simulating pigeons homing 
%%####behaviors. Homing pigeons can easily find their homes by using 
%%####three homing tools: magnetic field, sun and landmarks. In this  
%%#### newly invented algorithm, map and compass operator model is 
%%#### presented based on magnetic field and sun, while landmark operator 
%%#### model is presented based on landmarks. For some tough functions,  
%%#### it can quickly find the optimum, and it performs powerfully. For the 
%%#### most important reason, it combines some advantages of algorithms 
%%#### such as particle swarm optimization and artificial fish school algorithm. 
%***************initialization******************* 
clear all
clc 
format long;
T1=100;     %Global search algebra 
T2=100;     %Local search algebra 
pigeonnum=30;    %number  
D = 30;     % dimensionality 
R=0.3;     %parameters of magnetic field  
bound=[-40,40];    %hunting zone
tol = 1e-7;
f1=@fitness;
%**************initialization of the individual pigeon************ 
for i=1:pigeonnum                                                           %时间复杂度O(pigeonum*D*2)
    for j=1:D 
        x(i,j)=bound(1)+rand*(bound(2)-bound(1)); 
        v(i,j)=rand; 
    end 
end 
%**************calculate the fitness of pigeon*********** 
for i=1:pigeonnum                                                           %时间复杂度O(pigeonum*2)
    p(i)=f1(x(i,:),D); 
    p_best(i,:)=x(i,:); 
end 
%**************find the optimal pigeons******************** 
 
g_best=x(1,:); 
for i=2:pigeonnum                                                           %时间复杂度O(pigeonum-1)
    if f1(g_best,D)>f1(x(i,:),D) 
        g_best=x(i,:); 
    end 
end 
%************  magnetic compass and solar operator******************** 
for t=1:T1                                                                  %时间复杂度O(T1*(pigeonum*(2D+5))+1)
    for i=1:pigeonnum                                                       %时间复杂度O(pigeonum*(2D+5))                       
        v(i,:)=v(i,:)+rand*(p_best(i,:)-x(i,:)); 
        x(i,:)=x(i,:)*(1-exp(-R*t))+v(i,:);   %check whether beyond the searching space 
        for j=1:D                                    % magnetic field and solar operator 
            if abs(i-1)<=eps 
                if x(i,j)<bound(1)||x(i,j)>bound(2) 
                    x(i,j)=bound(1)+rand*(bound(2)-bound(1)); 
                    v(i,j)=rand; 
                end 
            else 
                if x(i,j)<bound(1)||x(i,j)>bound(2) 
                    x(i,j)=x(i-1,j);  
                    v(i,j)=v(i-1,j); 
                end     
            end 
        end 
        if f1(x(i,:),D)<p(i)                         %renewal individual fitness 
            p(i)=f1(x(i,:),D); 
            p_best(i,:)=x(i,:); 
        end 
        if p(i)<f1(g_best,D)                         %renewal global fitness 
            g_best=p_best(i,:); 
        end
    end
    result(t)=f1(g_best,D); 
end 
%*************???????********************** 
for t=1:T2                                                                  %时间复杂度O(T2*pigeonum*pigeonum)
    for i=1:pigeonnum-1                             %sort the pigeons       
        for j=i+1:pigeonnum 
            if f1(x(i,:),D)>f1(x(j,:),D) 
                temp_pigeon=x(i,:); 
                x(i,:)=x(j,:); 
                x(j,:)=temp_pigeon; 
            end 
        end 
    end 
    pigeonnum=ceil(pigeonnum/2);               %remove half of the pigeons according to the landmark 
    addpigeonnum=0;                        
    for i=1:pigeonnum 
        addpigeonnum=addpigeonnum+x(i,:);      
        p(i)=f1(x(i,:),D);                     %calculate fitness and location of the pigeon after sorting 
        p_best(i,:)=x(i,:); 
    end 
    pigeoncenter=ceil(addpigeonnum./pigeonnum);%calculate central position 
    for i=1:pigeonnum                                %local searching        
        for j=1:D                                    %check whether beyond the searching space 
            x(i,j) = x(i,j) + rand*(pigeoncenter(j)-x(i,j));
            while x(i,j)<bound(1)||x(i,j)>bound(2)
                x(i,j) = x(i,j) + rand*(pigeoncenter(j)-x(i,j));
            end
        end 
        if f1(x(i,:),D)<p(i)                         %renewal individual fitness 
            p(i)=f1(x(i,:),D); 
            p_best(i,:)=x(i,:); 
        end 
        if p(i)<f1(g_best,D)                         %renewal global fitness 
            g_best=p_best(i,:); 
        end 
    end 
    result(t+T1)=f1(g_best,D)
end 

figure                                               %graph 
for t=1:T1+T2-1 
    plot([t,t+1],[result(t),result(t+1)]); 
    hold on; 
end
xlabel('Iterative curve');
ylabel('Function value');
title(['The best value is ' num2str(min(result))])

end