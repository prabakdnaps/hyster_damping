nonl =  importdata('E:\Thesis\aorta\soft.txt','\t');
plot(nonl(:,1),nonl(:,2))
X=zeros(length(nonl),2);
X(:,2)=nonl(:,1);
Y=zeros(length(nonl),2);
Y(:,2)=nonl(:,2);
[val1,minidx]=min(X(:,2));
[val2,maxidx]=max(X(:,2));
val3=min(Y(:,2));
%finding the maximum and minimum values for splitting the curves in to
%two halves
index_sort=[minidx,maxidx];
index_sort=sort(index_sort);
x1=X(index_sort(1):index_sort(2),2);
y1=Y(index_sort(1):index_sort(2),2);
x2=X(1:index_sort(1),2);
x3=X(index_sort(2):end,2);
x2=[x2;x3];
y2=Y(1:index_sort(1),2);
y3=Y(index_sort(2):end,2);
y2=[y2;y3];
%bringing all values to positive side
x3=x1-val1;
y3=y1-val3;
x4=x2-val1;
y4=y2-val3;
[values,order]=sort(x4);
x4=x4(order,:);
y4=y4(order,:);
[values,order]=sort(x3);
x3=x3(order,:);
y3=y3(order,:);
%area of the two halves using trapizoidal method
area=abs(trapz(x3,y3)-trapz(x4,y4));
%find the middle line
low=[x1,y1];
high=[x2,y2];
[val,idx]=sort(low(:,1));
low=low(idx,:);
[val,idx]=sort(high(:,1));
high=high(idx,:);
high_mod=zeros(length(low),2);
min_length=min(length(low),length(high));
for i=1:min_length
    if low(i,1)==high(i,1)
        high_mod(i,:)=high(i,:);
    else
        idx=find(high(:,1)<=low(i,1));
        idx=idx(end);
        first_point=high(idx,:);
        if idx+1>length(high)
            second_point=high(idx,:);
        else
            second_point=high(idx+1,:);
        end
        high_mod(i,1)=low(i,1);
        m=((second_point(2)-first_point(2))/(second_point(1)-first_point(1)));
        high_mod(i,2)=(m*low(i,1))+first_point(2)-(m*first_point(1));
    end
end
high_mod(length(low),:)=high(end,:);
mid_line=[low(:,1),(low(:,2)+high_mod(:,2))/2];
k1=mid_line(:,1)\mid_line(:,2);
data=k1.*mid_line(:,1);
plot(mid_line(:,1),data,'-');
xlim([-0.1,0.1]);
ylim([-0.1,0.1]);
hold on;
plot(low(:,1),low(:,2),'+');
plot(high(:,1),high(:,2),'+');
plot(mid_line(:,1),mid_line(:,2),'+');
pause(5);
close all;