load('x1.mat');
load('x2.mat');
load('y1.mat');
load('y2.mat');
low=[x1,y1];
high=[x2,y2];
[val,idx]=sort(low(:,1));
low=low(idx,:);
[val,idx]=sort(high(:,1));
high=high(idx,:);
high_mod=zeros(length(high),2);
for i=1:length(high)
    if low(i,1)==high(i,1)
        high_mod(i,:)=high(i,:);
    else
        idx=find(high(:,1)<low(i,1));
        idx=idx(end);
        first_point=high(idx,:);
        second_point=high(idx+1,:);
        high_mod(i,1)=low(i,1);
        m=((second_point(2)-first_point(2))/(second_point(1)-first_point(1)));
        high_mod(i,2)=(m*low(i,1))+first_point(2)-(m*first_point(1));
    end
end
high_mod(length(low),:)=high(end,:);
mid_line=[low(:,1),(low(:,2)+high_mod(:,2))/2];
data=mid_line(:,1)\mid_line(:,2).*mid_line(:,1);
plot(mid_line(:,1),data,'-');
hold on;
plot(low(:,1),low(:,2),'+')
plot(high(:,1),high(:,2),'+')
plot(mid_line(:,1),mid_line(:,2),'+')