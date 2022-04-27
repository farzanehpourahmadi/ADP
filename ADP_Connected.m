tic
% clear
% %load('data')
% ro(1:5000,1:300)=20;
% phi1=0.6874;
% teta1=0.9234;
% teta24=0.8502;
% teta168=0.9665;
% for n=1:1000
% ep=normrnd(0,0.2369,[1000,1]);
% % ep=normrnd(0,0,[1000,1]);
% 
% for t=200:300
%     ro(n,t)=ro(n,t-1)+phi1*ro(n,t-1)-phi1*ro(n,t-2)+ro(n,t-24)-ro(n,t-25)-phi1*ro(n,t-25)+phi1*ro(n,t-26)...
%     +ro(n,t-168)-ro(n,t-169)-phi1*ro(n,t-169)+phi1*ro(n,t-170)-ro(n,t-192)+ro(n,t-193)+phi1*ro(n,t-193)-phi1*ro(n,t-194)...
%     +ep(t)-teta1*ep(t-1)-teta24*ep(t-24)+teta1*teta24*ep(t-25)-teta168*ep(t-168)+teta1*teta168*ep(t-169)+teta24*teta168*ep(t-192)-teta1*teta24*teta168*ep(t-193);
% 
%     if ro(n,t)<0
%        ro(n,t)=3;
%     end
% end
% end
% 
% % for n=2:1000
% %     ro(n,:)=ro(1,:);
% % end  
% 
% phi1_v=[0.9899,0.9775];
% teta1_v=[1.3156,1.4442];
% teta2_v=[-0.3504,-0.5509];
% teta41_v=[0.8424,0.8304];
% 
% 
% for n=1:1000
%     v{n}(1:1000,1:6)=50;
% eps=mvnrnd([0 0 0 0 0 0],[0.4288,0,0,0,0,0;0,0.0271,0,0,0,0;0,0,0.4288,0,0,0;0,0,0,0.0271,0,0;0,0,0,0,0.42889,0;0,0,0,0,0,0.02713],1000);
% % eps=mvnrnd([0 0 0 0 0 0],zeros(6,6),1000);
% 
% for t=50:300
%    v{n}(t,1)=v{n}(t-1,1)+phi1_v(1)*v{n}(t-1,1)-phi1_v(1)*v{n}(t-2,1)+eps(t,1)-teta1_v(1)*eps(t-1,1)-teta2_v(1)*eps(t-2,1)-teta41_v(1)*eps(t-41,1)+teta1_v(1)*teta41_v(1)*eps(t-42,1)+teta2_v(1)*teta41_v(1)*eps(t-43,1); 
%    v{n}(t,2)=v{n}(t-1,2)+phi1_v(2)*v{n}(t-1,2)-phi1_v(2)*v{n}(t-2,2)+eps(t,2)-teta1_v(2)*eps(t-1,2)-teta2_v(2)*eps(t-2,2)-teta41_v(2)*eps(t-41,2)+teta1_v(2)*teta41_v(2)*eps(t-42,2)+teta2_v(2)*teta41_v(2)*eps(t-43,2);
%    v{n}(t,3)=v{n}(t-1,3)+phi1_v(1)*v{n}(t-1,3)-phi1_v(1)*v{n}(t-2,3)+eps(t,3)-teta1_v(1)*eps(t-1,3)-teta2_v(1)*eps(t-2,3)-teta41_v(1)*eps(t-41,3)+teta1_v(1)*teta41_v(2)*eps(t-42,3)+teta2_v(1)*teta41_v(1)*eps(t-43,3);
%    v{n}(t,4)=v{n}(t-1,4)+phi1_v(2)*v{n}(t-1,4)-phi1_v(2)*v{n}(t-2,4)+eps(t,4)-teta1_v(2)*eps(t-1,4)-teta2_v(2)*eps(t-2,4)-teta41_v(2)*eps(t-41,4)+teta1_v(2)*teta41_v(2)*eps(t-42,4)+teta2_v(2)*teta41_v(2)*eps(t-43,4);
%    v{n}(t,5)=v{n}(t-1,5)+phi1_v(1)*v{n}(t-1,5)-phi1_v(1)*v{n}(t-2,5)+eps(t,5)-teta1_v(1)*eps(t-1,5)-teta2_v(1)*eps(t-2,5)-teta41_v(1)*eps(t-41,5)+teta1_v(1)*teta41_v(2)*eps(t-42,5)+teta2_v(1)*teta41_v(1)*eps(t-43,5);
%    v{n}(t,6)=v{n}(t-1,6)+phi1_v(2)*v{n}(t-1,6)-phi1_v(2)*v{n}(t-2,6)+eps(t,6)-teta1_v(2)*eps(t-1,6)-teta2_v(2)*eps(t-2,6)-teta41_v(2)*eps(t-41,6)+teta1_v(2)*teta41_v(2)*eps(t-42,6)+teta2_v(2)*teta41_v(2)*eps(t-43,6);
%     if v{n}(t,1)<0
%        v{n}(t,1)=0;
%     end
%     if v{n}(t,2)<0
%        v{n}(t,2)=0;
%     end
% end
% end
% % 
% % for n=2:1000
% %     v{n}=v{1};
% % end


R1=zeros(6,10);R1(1,1)=1;R1(2,2)=1;R1(3,3)=1;R1(4,4:6)=1;R1(4,8)=1;R1(5,7)=1;R1(5,9)=1;R1(6,10)=1;
R2=zeros(6,10);R2(1,6)=1;R2(2,7)=1;R2(3,8)=1;R2(4,1)=1;R2(4,3)=1;R2(4,9:10)=1;R2(5,2)=1;R2(5,4)=1;R2(6,5)=1;

n_i=400; %number of iterations
n_t=48;
n_j=6; %number of reseviors
n_l=10; %number of lines

e=[0.1 0.04 0.03 0.1 0.03 -1.1 -1.04 -1.03 -1.1 -1.03];
lmax=[6.507;0.51;0.628;10.73;2.853;0.514]*10;
lmin=0.1*lmax;
linit=0.11*lmax;
pdmax=[0.2397;0.114;0.0228;0.3022;0.0239;0.0114]*10;
pdmin=0*pdmax;
pmax=[0.2522;0.2522;0.0228;0.2522;0.2522;0.2522;0.2522;0.0228;0.2522;0.2522];

for t=1:n_t
a{1,t}=[ro(2,200+t)*e(1),ro(2,200+t)*e(2),ro(2,200+t)*e(3),ro(2,200+t)*e(4),ro(2,200+t)*e(5),ro(2,200+t)*e(5)];
b{1,t}=[t*ro(2,200+t)*e(1),t*ro(2,200+t)*e(2),t*ro(2,200+t)*e(3),t*ro(2,200+t)*e(4),t*ro(2,200+t)*e(5),t*ro(2,200+t)*e(5)];
% a{1,t}=[(sum(ro(2:n_i,200+t))/(n_i-1))*e(1),(sum(ro(2:n_i,200+t))/(n_i-1))*e(2)];
end

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     00;
for n=501:600
    n
     for t=2:n_t
    pd=sdpvar(n_j,1,'full');
    pc=sdpvar(n_j,1,'full');
    p=sdpvar(n_l,1,'full');
    l=sdpvar(n_j,1,'full');
    lp=sdpvar(n_j,1,'full');
    Constraints=[];
    for j=1:n_j     
%         Constraints=[Constraints,l(j)<=lmax(j)];
        Constraints=[Constraints,l(j)>=lmin(j)];
%         Constraints=[Constraints,l(j)<=lmax(j)];
        Constraints=[Constraints,pd(j)<=pdmax(j)];
        Constraints=[Constraints,pd(j)>=pdmin(j)];
        Constraints=[Constraints,pc(j)<=pdmax(j)];
        Constraints=[Constraints,pc(j)>=pdmin(j)];
        Constraints=[Constraints,pd(j)==R1(j)*p(j)];
        Constraints=[Constraints,pc(j)==R2(j)*p(j)];
    end
    for j=1:n_l
        Constraints=[Constraints,p(j)<=pmax(j)];
        Constraints=[Constraints,p(j)>=0];
    end
        Constraints=[Constraints,l==linit-pd+0.6*pc+(v{n}(t+50,:)/500)'];
        Constraints=[Constraints,lp==linit-pd+0.6*pc];

    Objective=ro(n,t+200)*e*p+a{n-1,t}*lp+b{n-1,t}*(v{n}(t+50,:)/500)';

    ops = sdpsettings('solver','gurobi','verbose',0);
    sol=optimize(Constraints,-Objective,ops);
c=value(-Objective);
l1{n,t}=value(l);
p1{n,t}=value(p);
VE{n,t}=c;

for jj=1:n_j
    uni=zeros(n_j,1);
    uni(jj)=1;
    pd=sdpvar(n_j,1,'full');
    pc=sdpvar(n_j,1,'full');
    p=sdpvar(n_l,1,'full');
    l=sdpvar(n_j,1,'full');
    lp=sdpvar(n_j,1,'full');
    Constraints=[];
    for j=1:n_j     

        Constraints=[Constraints,l(j)>=lmin(j)];
%         Constraints=[Constraints,l(j)<=lmax(j)];
        Constraints=[Constraints,pd(j)<=pdmax(j)];
        Constraints=[Constraints,pd(j)>=pdmin(j)];
        Constraints=[Constraints,pc(j)<=pdmax(j)];
        Constraints=[Constraints,pc(j)>=pdmin(j)];
        Constraints=[Constraints,pd(j)==R1(j)*p(j)];
        Constraints=[Constraints,pc(j)==R2(j)*p(j)];
        Constraints=[Constraints,l==linit+uni-pd+0.6*pc+(v{n}(t+50,:)/500)'];
        Constraints=[Constraints,lp==linit+uni-pd+0.6*pc];
    end
    for j=1:n_l
        Constraints=[Constraints,p(j)<=pmax(j)];
        Constraints=[Constraints,p(j)>=0];
    end

        Objective=ro(n,t+200)*e*p+a{n-1,t}*lp+b{n-1,t}*(v{n}(t+50,:)/500)';
%         Objective=ro(n,t+200)*p(j)+Objective;

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
            pd=sdpvar(n_j,1,'full');
            pc=sdpvar(n_j,1,'full');
            p=sdpvar(n_l,1,'full');
            l=sdpvar(n_j,1,'full');
            lp=sdpvar(n_j,1,'full');
            Constraints=[];
    for j=1:n_j     

        Constraints=[Constraints,l(j)>=lmin(j)];
%         Constraints=[Constraints,l(j)<=lmax(j)];
        Constraints=[Constraints,pd(j)<=pdmax(j)];
        Constraints=[Constraints,pd(j)>=pdmin(j)];
        Constraints=[Constraints,pc(j)<=pdmax(j)];
        Constraints=[Constraints,pc(j)>=pdmin(j)];
        Constraints=[Constraints,pd(j)==R1(j)*p(j)];
        Constraints=[Constraints,pc(j)==R2(j)*p(j)];
        Constraints=[Constraints,l==linit+uni-pd+0.6*pc+(v{n}(t+50,:)/500)'];
        Constraints=[Constraints,lp==linit+uni-pd+0.6*pc];
    end
    for j=1:n_l
        Constraints=[Constraints,p(j)<=pmax(j)];
        Constraints=[Constraints,p(j)>=0];
    end
            
                Objective=ro(n,t+200)*e*p+a{n-1,t}*lp+b{n-1,t}*((v{n}(t+50,:)/500)'+uni);
%               Objective=ro(n,t+200)*p(j)+Objective;

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
a{n,n_t}=[ro(n,200+n_t+1)*e(1),ro(n,200+n_t+1)*e(2),ro(n,200+n_t+1)*e(3),ro(n,200+n_t+1)*e(4),ro(n,200+n_t+1)*e(5),ro(n,200+n_t+1)*e(5)];
b{n,n_t}=[n_t*ro(n,200+n_t+1)*e(1),n_t*ro(n,200+n_t+1)*e(2),n_t*ro(n,200+n_t+1)*e(3),n_t*ro(n,200+n_t+1)*e(4),n_t*ro(n,200+n_t+1)*e(5),n_t*ro(n,200+n_t+1)*e(5)];
n
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
    l13(n,t)=l1{n,t}(3);
    l14(n,t)=l1{n,t}(4);
    l15(n,t)=l1{n,t}(5);
    l16(n,t)=l1{n,t}(6);
    end
end


for n=2:n_i
for t=2:n_t
CC1(n,t)=VE{n,t};
end
end
