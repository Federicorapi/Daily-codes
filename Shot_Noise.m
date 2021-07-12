clc
clear all
close all
%% Laser
Laser = zeros(2000,2,1000);
Alexa = zeros(2000,2,1000);
for i  = 1:9
    T = readtable(strcat('laser/laser000',num2str(i),'.asc'),'filetype','text');
    Laser(:,:,i) = table2array(T);
end
for i  = 10:99
    T = readtable(strcat('laser/laser00',num2str(i),'.asc'),'filetype','text');
    Laser(:,:,i) = table2array(T);
end

for i  = 100:999
    T = readtable(strcat('laser/laser0',num2str(i),'.asc'),'filetype','text');
    Laser(:,:,i) = table2array(T);
end

i = 1000;
    T = readtable(strcat('laser/laser',num2str(i),'.asc'),'filetype','text');
    Laser(:,:,i) = table2array(T);
  
for i = 1:1000
    background = mean(Laser(1:926,2,i));
Z_l(i) = trapz(Laser(1000:1007,2,i)-background);
end

[mean_laser,std_laser] = normfit(Z_l);
shot_noise_laser = sqrt(mean_laser);


%% Alexa
for i  = 1:9
    T = readtable(strcat('alexa/alexa000',num2str(i),'.asc'),'filetype','text');
    Alexa(:,:,i) = table2array(T);
end
for i  = 10:99
    T = readtable(strcat('alexa/alexa00',num2str(i),'.asc'),'filetype','text');
    Alexa(:,:,i) = table2array(T);
end

for i  = 100:999
    T = readtable(strcat('alexa/alexa0',num2str(i),'.asc'),'filetype','text');
    Alexa(:,:,i) = table2array(T);
end

i = 1000;
    T = readtable(strcat('alexa/alexa',num2str(i),'.asc'),'filetype','text');
    Alexa(:,:,i) = table2array(T);
  
for i = 1:1000
    background = mean(Alexa(1:599,2,i));
Z_a(i) = trapz(Alexa(1100:1107,2,i)-background);
end

[mean_alexa,std_alexa] = normfit(Z_a);
shot_noise_alexa = sqrt(mean_alexa);



figure, hold on
plot(Z_a)
plot(Z_l)


%% Kinetics image mode
for i  = 1:9
    T = readtable(strcat('luminescence/luminescence000',num2str(i),'.asc'),'filetype','text');
    image(:,:,i) = table2array(T(:,2:81));
end
for i  = 10:99
    T = readtable(strcat('luminescence/luminescence00',num2str(i),'.asc'),'filetype','text');
    image(:,:,i) = table2array(T(:,2:81));
end

for i  = 100:300
    T = readtable(strcat('luminescence/luminescence0',num2str(i),'.asc'),'filetype','text');
    image(:,:,i) = table2array(T(:,2:81));
end

%% Work on the image

for i = 1:2000
    for j = 1:80
        mean_image(i,j) = mean(image(i,j,1:300));
    end
end


image2plot = image-mean_image; 

figure, imagesc(image2plot(900:1100,:,1)');
colorbar
gif('image_gif.gif','DelayTime',100,'overwrite',true)

for i = 2:300
 imagesc(image2plot(900:1100,:,i)');
 colorbar
 gif
 
end

%% Long kinetics
for i  = 1:9
    T = readtable(strcat('long kinetics/long kinetics000',num2str(i),'.asc'),'filetype','text');
    spectro(:,:,i) = table2array(T);
end
for i  = 10:99
    T = readtable(strcat('long kinetics/long kinetics00',num2str(i),'.asc'),'filetype','text');
    spectro(:,:,i) = table2array(T);
end

for i  = 100:999
    T = readtable(strcat('long kinetics/long kinetics0',num2str(i),'.asc'),'filetype','text');
    spectro(:,:,i) = table2array(T);
end

for i  = 1000:4000
    T = readtable(strcat('long kinetics/long kinetics',num2str(i),'.asc'),'filetype','text');
    spectro(:,:,i) = table2array(T);
end

for i = 1:4000
    background = mean(spectro(1:571,2,i));
    Z(i) = trapz(spectro(735:1151,2,i)-background);
end
figure, plot(Z)
grid on