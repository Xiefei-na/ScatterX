folder_path = uigetdir();


file_list1 = dir(fullfile(folder_path, '*.tif'));  % 获取所有Iochamber格式的文件
P = cell(length(file_list1), 1);  % 创建一个单元数组来存储文件名
for i = 1:length(file_list1)
P{i} = [folder_path,'\',file_list1(i).name];  % 将文件名存储在P中
end



%input background image and get baemcenter
[fileName, filePath] = uigetfile('*.tif');
p = double(imread([filePath,fileName]));
p(p<1)=NaN;p(p>100000)=NaN;
[cenx,ceny] = getBeamCenter(p,6);
imagesc(p,[0 500]);hold on;scatter(cenx,ceny,'r','filled');




% cenx=gpuArray(401);    ceny=291;
phi1=136;    phi2=225;
rmax=1000;   sampletodetector=1139;    pixelsize=0.172;
lambda=0.154;
[qi,Ii] = arrayfun(@(i) cake2qintensity(double(imread(P{i})),cenx,ceny,phi1,phi2,rmax,lambda,sampletodetector,pixelsize),(1:i),'UniformOutput',false) ;



q1=0.3;
q2=2.2;
[qusefulli,Iusefulli] = arrayfun(@(i) cutuseless(qi{i},Ii{i},q1,q2),(1:i),'UniformOutput',false) ;

Porodlimit = 4;
Guinierlimit = 0.4;
Exlimit = 10;
[qfulli,Ifulli] = arrayfun(@(i) extSAXS(qusefulli{i}',Iusefulli{i}',Porodlimit,Guinierlimit,Exlimit),(1:i),'UniformOutput',false) ;

tic
[r,Vr] = arrayfun(@(i) chahine2pdf(qfulli{i}',Ifulli{i}'),(1:i),'UniformOutput',false) ;
toc


[fitresult, gof] = arrayfun(@(i) createFit(qusefulli{i},Iusefulli{i}),(1:i),'UniformOutput',false) ;
