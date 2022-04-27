tic
% clear
% %load('data')
% ro(1:1000)=174;
% phi1=0.6874;
% teta1=0.9234;
% teta24=0.8502;
% teta168=0.9665;
% % ep=normrnd(0,2.4369,[1000,1]);
% ep=normrnd(0,0,[1000,1]);
% 
% for t=200:1000
%     ro(t)=ro(t-1)+phi1*ro(t-1)-phi1*ro(t-2)+ro(t-24)-ro(t-25)-phi1*ro(t-25)+phi1*ro(t-26)...
%     +ro(t-168)-ro(t-169)-phi1*ro(t-169)+phi1*ro(t-170)-ro(t-192)+ro(t-193)+phi1*ro(t-193)-phi1*ro(t-194)...
%     +ep(t)-teta1*ep(t-1)-teta24*ep(t-24)+teta1*teta24*ep(t-25)-teta168*ep(t-168)+teta1*teta168*ep(t-169)+teta24*teta168*ep(t-192)-teta1*teta24*teta168*ep(t-193);
% end  
% 
% phi1_v=[0.9899,0.9775];
% teta1_v=[1.3156,1.4442];
% teta2_v=[-0.3504,-0.5509];
% teta41_v=[0.8424,0.8304];
% v(1:1000,1:2)=20;
% % eps=mvnrnd([0 0],[8.577924,0.090712;0.090712,0.542625],1000);
% eps=mvnrnd([0 0],[0 0;0 0],1000);
% 
% for t=50:1000
%    v(t,1)=v(t-1,1)+phi1_v(1)*v(t-1,1)-phi1_v(1)*v(t-2,1)+eps(t,1)-teta1_v(1)*eps(t-1,1)-teta2_v(1)*eps(t-2,1)-teta41_v(1)*eps(t-41,1)+teta1_v(1)*teta41_v(1)*eps(t-42,1)+teta2_v(1)*teta41_v(1)*eps(t-43,1); 
%    v(t,2)=v(t-1,2)+phi1_v(2)*v(t-1,2)-phi1_v(2)*v(t-2,2)+eps(t,2)-teta1_v(2)*eps(t-1,2)-teta2_v(2)*eps(t-2,2)-teta41_v(2)*eps(t-41,2)+teta1_v(2)*teta41_v(2)*eps(t-42,2)+teta2_v(2)*teta41_v(2)*eps(t-43,2); 
% end


n_i=100; %number of iterations
n_t=48; 
% e=[0.1101 0.5051];
% lmax=[11300;2500];
% lmin=0.1*lmax;
% linit=0.9*lmax;
% pmax=[5.796;6.336];
% pmin=0.1*pmax;

% e=[0.1101 0.5051];
e=[0.5051 0.5051];
lmax=0.1*[11300;10000];
lmin=0.1*lmax;
linit=0.11*lmax;
pmax=10*[5.796;12.136];
pmin=0*pmax;
bc=[2;1];

for t=1:n_t
a{1,t}=[ro(2,200+t)*e(1),ro(2,200+t)*e(2)];
b{1,t}=[2*n_t*ro(2,200+t)*e(1),n_t*ro(2,200+t)*e(2)];
% a{1,t}=[(sum(ro(2:n_i,200+t))/(n_i-1))*e(1),(sum(ro(2:n_i,200+t))/(n_i-1))*e(2)];

end
n_j=2; %number of reseviors
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     00;
for n=2:n_i
    n
     for t=2:n_t
    p=sdpvar(n_j,1,'full');
    l=sdpvar(n_j,1,'full');
    lp=sdpvar(n_j,1,'full');
    Constraints=[];
    for j=1:n_j     
%         Constraints=[Constraints,l(j)<=lmax(j)];
        Constraints=[Constraints,l(j)>=lmin(j)]; 
        Constraints=[Constraints,p(j)<=pmax(j)];
        Constraints=[Constraints,p(j)>=pmin(j)];
        Constraints=[Constraints,l==linit-[1 0;-1 1]*p+v{n}(t+50,:)'];
        Constraints=[Constraints,lp==linit-[1 0;-1 1]*p];
    end
    Objective=0;
    for j=1:n_j
            Objective=ro(n,t+200)*e(j)*p(j)+a{n-1,t}(j)*lp(j)+b{n-1,t}(j)*v{n}(t+50,j)+Objective;
%           Objective=ro(n,t+200)*p(j)+Objective;
    end
    ops = sdpsettings('solver','gurobi','verbose',0);
    sol=optimize(Constraints,-Objective,ops);
c=value(-Objective);
l1{n,t}=value(l);
p1{n,t}=value(p);
VE{n,t}=c;

for jj=1:n_j
    uni=zeros(n_j,1);
    uni(jj)=1;
    p=sdpvar(n_j,1,'full');
    l=sdpvar(n_j,1,'full');
    lp=sdpvar(n_j,1,'full');
    Constraints=[];
    for j=1:n_j     
%     Constraints=[Constraints,l(j)<=lmax(j)];
    Constraints=[Constraints,l(j)>=lmin(j)]; 
    Constraints=[Constraints,p(j)<=pmax(j)];
    Constraints=[Constraints,p(j)>=pmin(j)];
    Constraints=[Constraints,l==linit+uni-[1 0;-1 1]*p+v{n}(t+50,:)'];
    Constraints=[Constraints,lp==linit+uni-[1 0;-1 1]*p];

    end
    Objective=0;
    for j=1:n_j
        Objective=ro(n,t+200)*p(j)*e(j)+a{n-1,t}(j)*lp(j)+b{n-1,t}(j)*v{n}(t+50,j)+Objective;
%         Objective=ro(n,t+200)*p(j)+Objective;
    end
    ops = sdpsettings('solver','gurobi','verbose',0);
    sol=optimize(Constraints,-Objective,ops);
    cc(jj)=value(-Objective);
l2{n,t}=value(l);
p2{n,t}=value(p);
a{n,t-1}(jj)=0.5*a{n-1,t-1}(jj)+0.5*(-cc(jj)+c);
end

        for jj=1:n_j
            uni=zeros(n_j,1);
            uni(jj)=1;
            p=sdpvar(n_j,1,'full');
            l=sdpvar(n_j,1,'full');
            lp=sdpvar(n_j,1,'full');
            Constraints=[];
            for j=1:n_j     
%               Constraints=[Constraints,l(j)<=lmax(j)];
                Constraints=[Constraints,l(j)>=lmin(j)]; 
                Constraints=[Constraints,p(j)<=pmax(j)];
                Constraints=[Constraints,p(j)>=pmin(j)];
                Constraints=[Constraints,l==linit+uni-[1 0;-1 1]*p+v{n}(t+50,:)'];
                Constraints=[Constraints,lp==linit-[1 0;-1 1]*p];

            end
            Objective=0;
            for j=1:n_j
                Objective=ro(n,t+200)*p(j)*e(j)+a{n-1,t}(j)*lp(j)+b{n-1,t}(j)*(v{n}(t+50,j)+uni(j))+Objective;
%               Objective=ro(n,t+200)*p(j)+Objective;
            end
            ops = sdpsettings('solver','gurobi','verbose',0);
            sol=optimize(Constraints,-Objective,ops);
            cc1(jj)=value(-Objective);
            l2{n,t}=value(l);
            p2{n,t}=value(p);
            b{n,t-1}(jj)=0.5*b{n-1,t-1}(jj)+0.5*(-cc1(jj)+c);
        end


linit=l1{n,t};
     end
linit=0.11*lmax;
% a{n,n_t}=a{n,n_t-1};  
% a{n,n_t}=[(sum(ro(n,200+n_t+1:200+2*n_t))/(200+2*n_t-(200+n_t+1)))*e(1),(sum(ro(n,200+n_t+1:200+2*n_t))/(200+2*n_t-(200+n_t+1)))*e(2)];
a{n,n_t}=[ro(n,200+n_t+1)*e(1),ro(n,200+n_t+1)*e(2)];
b{n,n_t}=[2*n_t*ro(n,200+n_t+1)*e(1),n_t*ro(n,200+n_t+1)*e(2)];
end
 time=toc
 
 %% result
 for n=2:n_i
    for t=2:n_t
    aa1(n,t)=a{n,t}(1);
    aa2(n,t)=a{n,t}(2);
    end
 end

 for n=2:n_i
    for t=2:n_t
    bb1(n,t)=b{n,t}(1);
    bb2(n,t)=b{n,t}(2);
    end
end

for n=2:n_i
    for t=2:n_t
    p11(n,t)=p1{n,t}(1);
    p12(n,t)=p1{n,t}(2);
    end
end

for n=2:n_i
    for t=2:n_t
    l11(n,t)=l1{n,t}(1);
    l12(n,t)=l1{n,t}(2);
    end
end

for n=2:n_i
    for t=2:n_t
    Value(n,t)=ro(n,t+200)*e(1)*p11(n,t)+ro(n,t+200)*e(2)*p12(n,t);
    end
end

for n=2:n_i
for t=2:n_t
CC1(n,t)=VE{n,t};
end
end

v_c(1:2,1:48)=0;
for n=3:500
for t=2:47
v_c(n,t)=-VE{n,t}+v_c(n-1,t+1)-a{n-1,t}*(l1{n-1,t}-(v{n}(t+50,:))')-b{n-1,t}*(v{n}(t+50,:))';
end
end
