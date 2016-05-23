
clear all;

%%%%%%%%%%%%%%% model selection %%%%%%%%%%%%%%%
method=5;
% 1 = random (no learning)
% 2 = CPM
% 3 = Arthur's model
% 4 = Roth & Erev's model
% 5 = apiration

%parameters in Arthur's model
C=1.1;
p=0.5;

%parameters in Roth & Erev model
lambda=0.90;

%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%

w=100;%width
h=100;%height
N=8;%nr of skaters
k=60;%minimun difference in directions
delta=25;%speed
r=10;%collision distance

R1=1;%positive reward
R2=-0.8;%negative reward

%%%%%%%%%%%%%% initialization %%%%%%%%%%%%%%%%%%

n=360/k;%nr of actions

A=0:k:(n-1)*k;%set of directions/actions

R=zeros(N,1);%reward vector for all skaters
theta=zeros(N,n);
q(:,:,1)=ones(N,n)/n;%q at time=1

% positions
xpos=w*rand(N,1);
ypos=h*rand(N,1);


%%%%%%%%%%%%%%%%%%%%% start game %%%%%%%%%%%%%%%%%%%%%%%%%

for t=1:2000%2000 rounds
    for i=1:N% N skaters
        
        qtemp=q(i,:,t);
        if all(qtemp==0)
            qtemp(:)=1/n;
            a(i,t)=randsrc(1,1,[A;qtemp]);
        elseif all(qtemp>=0)
            qtemp=qtemp/sum(qtemp);
            a(i,t)=randsrc(1,1,[A;qtemp]);
        else%which means negative elements in q(i,:,t)
            qtemp=qtemp-min(qtemp);
%                 qtemp(qtemp(:)<0)=0;
            qtemp=qtemp/sum(qtemp);
            a(i,t)=randsrc(1,1,[A;qtemp]);
        end
        
        %hypothetical position of current skater
        hxpos=mod(xpos(i)+delta*cos(a(i,t)*2*pi/360),w);
        hypos=mod(ypos(i)+delta*sin(a(i,t)*2*pi/360),h);
        
        %distance
        xdis=min(abs(xpos-hxpos), abs(w-(xpos-hxpos)));%(vector)distance in x-axis
        ydis=min(abs(ypos-hypos), abs(h-(ypos-hypos)));%distance in y-axis
        dis=sqrt(xdis.^2+ydis.^2);
        
        if all(dis(:)>r*ones(N,1))
            xpos(i)=hxpos;
            ypos(i)=hypos;
            R(i,t)=R1;
        else
            R(i,t)=R2;
        end
        
        switch method
            case 1%random
                q(i,:,t+1)=ones(1,n)/n;
            case 2%CPM
                %calculate theta and q
                theta(i,:)=theta(i,:)+(a(i,t)*ones(1,n)==A)*R(i,t);%CPM, only care about time=t
                if t<50
                    q(i,:,t+1)=ones(1,n)/n;
                else
                    q(i,:,t+1)=theta(i,:)/sum(theta(i,:));
                end
            case 3%Arthur
                theta(i,:)=theta(i,:)+(a(i,t)*ones(1,n)==A)*R(i,t);
                if t<50
                    q(i,:,t+1)=ones(1,n)/n;
                else
                    dq=(R(i,t)/(R(i,t)+C*t^p))*((a(i,t)*ones(1,n)==A)-q(i,:,t));
                    q(i,:,t+1)=q(i,:,t)+dq;
                end
                
            case 4%R&E
                theta(i,:)=lambda*theta(i,:)+(a(i,t)*ones(1,n)==A)*R(i,t);
                if t<50
                    v(t)=sum(R(i,:));
                    q(i,:,t+1)=ones(1,n)/n;
                else
                    v(t)=lambda*v(t-1)+R(i,t);
                    dq=(R(i,t)/v(t))*((a(i,t)*ones(1,n)==A)-q(i,:,t));
                    q(i,:,t+1)=q(i,:,t)+dq;
                end  
            case 5%aspiration, average past payoffs
                theta(i,:)=theta(i,:)+(a(i,t)*ones(1,n)==A)*R(i,t);
                if t<50
                    q(i,:,t+1)=ones(1,n)/n;
                else
                    asp(t)=sum(theta(i,:))/t;
                    dq=(R(i,t)-asp(t))*((a(i,t)*ones(1,n)==A)-q(i,:,t));
                    q(i,:,t+1)=q(i,:,t)+dq;
                end
                
        end
    end

    for k = 1:n
        meanR(t,k)= sum(sum(R.*(a==A(k)*ones(size(R)))))/sum(sum(a==A(k)*ones(size(R))));
    end
end

%%%%%%%%%%%%%%%%%%% display results %%%%%%%%%%%%%%%%%%%%

%calculate final score
scoreM=R;
scoreM(scoreM(:,:)<0)=0;
scoreM(scoreM(:,:)>0)=1;
for t=1:2000
    scoreT(t)=sum(sum(scoreM(:,1:t)))/numel(scoreM(:,1:t));
end

figure;
plot(scoreT);
xlabel('timestep');
ylabel('score on the time');
title('score for all actions');
% axis tight;
ylim([0 1]);

score=sum(sum(scoreM))/numel(scoreM);

figure;
plot(meanR);
legend([num2str(A')]);
xlabel('timestep');
ylabel('average payoff');
title('average payoff for each action');
 axis tight;
%ylim([0 1]);

figure;
q1=permute(q,[2 3 1]);
plot(1:2001,q1(1,:,1),1:2001,q1(2,:,1),1:2001,q1(3,:,1),1:2001,q1(4,:,1),1:2001,q1(5,:,1),1:2001,q1(6,:,1));
legend([num2str(A')]);
xlabel('timestep');
ylabel('q');
title('change of q for only skater 1');
ylim([-0.1 1.1]);
% axis tight;


figure;
plot(a')
xlabel('timestep');
ylabel('action');
title('action selection upon time for all skaters');
axis tight;
ylim([-50 350]);