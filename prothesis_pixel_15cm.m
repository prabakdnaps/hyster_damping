%%extracting the data from the pictures
%after pressurization
ima=imread('E:\Thesis\aorta\IMG_9590.CR2');
ima = imcrop(ima,[759 1491 4114 1210]);
ima = imrotate(ima,90);
imshow(ima);
figure;
% a=ima(:,:,1)+ima(:,:,2)+ima(:,:,3);
a=ima(:,:,1);
save('a.mat','a');
%before pressurization
ima=imread('E:\Thesis\aorta\IMG_9588.CR2');
ima = imcrop(ima,[759 1491 4114 1210]);
ima = imrotate(ima,90);
imshow(ima);
b=ima(:,:,1);
save('b.mat','b');
%% creating the profile

load('a.mat');
profile=zeros(4115,4);
for i=1:4115
    test=find(a(i,:,1)>120);
    profile(i,1)=test(1);
    profile(i,2)=test(end);
end
plot(profile(:,1))
hold on;
plot(profile(:,2))

load('b.mat');
for i=1:4115
    test=find(b(i,:,1)>120);
    profile(i,3)=test(1);
    profile(i,4)=test(end);
end
plot(profile(:,3))
hold on;
plot(profile(:,4));
save('profile.mat','profile');
%% 
load('profile.mat');
xaxis=linspace(1,4115,4115);
after_left=[xaxis',profile(:,1)];
after_right=[xaxis',profile(:,2)];
before_left=[xaxis',profile(:,3)];
before_right=[xaxis',profile(:,4)];
I1=find(profile(:,1)>100);
after_left=after_left(I1,2);
plot(after_left);
hold on;
I2=find(profile(:,2)<1200);
after_right=after_right(I2,2);
plot(after_right);

I3=find(profile(:,3)>100);
before_left=before_left(I3,2);
plot(before_left);

I4=find(profile(:,4)<1200);
before_right=before_right(I4,2);
plot(before_right);
%% 
dia_af=profile(:,2)-profile(:,1);
dia_be=profile(:,4)-profile(:,3);
change_dia=dia_af-dia_be;
per_elon=change_dia./dia_be;
I=find(per_elon<0.5&per_elon>0);
final_elong=per_elon(I);
plot(final_elong);
save('final_elon.mat','final_elong');