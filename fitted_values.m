function [area,min_Disp,max_Disp,KK,kk] = fitted_values(X,Y)
%fitted_values returns the area under curve
%   Detailed explanation goes here
temp_idx=X(:,1)~=0;
X=X(temp_idx,:);
temp_idx=Y(:,1)~=0;
Y=Y(temp_idx,:);
%since the force is having a negative sign
Y=[Y(:,1),Y(:,2)];
%remove the ratio
X(:,2)=X(:,2);
coeff=zeros(1,8);
if length(X)>4
    foco1=fit(X(:,1),X(:,2),'fourier1');
    coeff(1:4)=coeffvalues(foco1);
    foco2=fit(Y(:,1),Y(:,2),'fourier1');
    coeff(5:8)=coeffvalues(foco2);
    mm=coeff(4).^2*(-coeff(2).*cos(coeff(4).*X(:,1))-coeff(3).*sin(coeff(4).*X(:,1)));
    cm=coeff(4)*(-coeff(2).*sin(coeff(4).*X(:,1))+coeff(3).*cos(coeff(4).*X(:,1)));
    km=coeff(1)+coeff(2).*cos(coeff(4).*X(:,1))+coeff(3).*sin(coeff(4).*X(:,1));
    fm=coeff(6).*cos(coeff(8).*X(:,1))+coeff(7).*sin(coeff(8).*X(:,1));
    %omega_n=103.6*2*pi;
    omega_n=1;
    %zeta=0.0388;
    k1m=mm./(omega_n).^2;
    k2m=(2*1/omega_n).*cm;
    k3m=km;
    A=[k2m,km.^2,km.^3];
    b=-((k1m+k3m).*1)+fm;
    kk=A\b;
end
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
    %pca analysis
    XX=[X(:,2),Y(:,2)];
    sigma=(XX*XX')*(1/length(XX));
    [U,S,V]=svd(sigma);
    Z=U(:,1:2)'*XX;
    XX_red=U(:,1:2)*Z;
%     subplot(1,2,1);
%     plot(XX(:,1),XX(:,2));
%     hold on
%     plot(XX_red(:,1),XX_red(:,2))
%     hold off
    %for identifying k3 &k2 ; the phase information is removed
%     [val1,idx1]=max(X(:,2));
%     [val2,idx2]=max(Y(:,2));
%     if idx1>idx2
%         temp1=idx1-idx2;
%         Y(:,2)=circshift(Y(:,2),temp1);
%     else
%         temp1=idx2-idx1;
%         Y(:,2)=circshift(Y(:,2),-temp1);
%     end
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
    x4=x4(order,:)+1;
    y4=y4(order,:)+1;
    [values,order]=sort(x3);
    x3=x3(order,:)+1;
    y3=y3(order,:)+1;
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
    for i=1:length(low)
        if i<length(high)
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
    XX=mid_line(:,1);
    OS=ones(length(XX),1);
    A=[OS,XX,XX.^2,XX.^3];
    KK=A\mid_line(:,2);
    data=A*KK;
    min_Disp=min(X(:,2));
    max_Disp=max(X(:,2));
%     subplot(1,2,2);
%     plot(mid_line(:,1),data,'-');
%     hold on;
%     plot(low(:,1),low(:,2),'+');
%     plot(high(:,1),high(:,2),'+');
%     plot(mid_line(:,1),mid_line(:,2),'+');
%     hold off;
%     pause(20);
%     close all;
else
    area=0;
    min_Disp=0;
    max_Disp=0;
    KK=0;
end
end

