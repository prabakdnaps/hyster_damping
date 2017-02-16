% clear all;
% disp=load('disp.mat');
% force=load('force.mat');
%removing the zeros and cleaning the data
X=disp_matrix(:,399:400);
Y=force_matrix(:,399:400);
temp_idx=X(:,1)~=0;
X=X(temp_idx,:);
temp_idx=Y(:,1)~=0;
Y=Y(temp_idx,:);
%since the force is having a negative sign
Y=[Y(:,1),-Y(:,2)];

%finding number of cycles
temp=circshift(X(:,2),1);
temp(1)=0;
temp_idx=(X(:,2).*temp)<0;
temp_idx1=find(temp_idx);
%if the data contains less than one cycle, area=0 otherwise area is
%calculated
if length(temp_idx1)>2
    X=X(temp_idx1(1):temp_idx1(3),:);
    Y=Y(temp_idx1(1):temp_idx1(3),:);
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
%     figure;
%     plot(x1,y1,x2,y2);
else
    area=0;
end