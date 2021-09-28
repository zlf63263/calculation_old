clear all; close all; clc;
Aoled = 10 * (1e-3)^2;
Apd = 12.96 * (1e-3)^2;
r= (59.5) * (1e-3);

prefix='angle';


% load spectrum ---------------------------------------
plotSpec =figure('units','normalized','outerposition',[0 0 1 1]);
F = dir('*.IRR');
C = {'b','r','g','y','m','k','c',[.2 .5 1], [.2 .6 .2]};
figure(plotSpec); subplot(211)
for ii = 1:length(F)
    filename = cellstr(F(ii).name);    
     mat=[];
   fid = fopen(F(ii).name);
    tline = fgets(fid);
    line_number = 1;
    
    while ischar(tline)
       line_number = line_number + 1;
       oline = strsplit(strtrim(tline));
       if ~isnan(str2double(oline(1))) % if the first entry in that line is a valid number, then u take that line and process and append to mat
           tmp_mat=zeros(1,length(oline));
           for j=1:length(oline)
               tmp_mat(j) = str2double(oline(j));
           end
           mat = cat(1,mat,tmp_mat);
       end
        tline = fgets(fid);
    end
    fclose(fid);
    Intensity_lamda_theta(:,ii) = mat(:,2);
    lamda = mat(:,1) * 1e-9;    
    
    ARAW = (trapz(lamda,Intensity_lamda_theta(:,ii)))^-1;
    rho_lamda_thetaRAW = zeros(length(Intensity_lamda_theta),1);
    rho_lamda_thetaRAW = Intensity_lamda_theta * ARAW;
    hold on;
    plot(lamda*1e9,rho_lamda_thetaRAW(:,ii),'color',C{ii},'LineWidth',3.5);
    
end
legend('0','10','20','30','40','50','60','70','80');
xlabel('lamda'); ylabel('normalized spectrum');

hold off;
% pad w/ zeros for rho(90) = 0; --------------------------------
Intensity_lamda_theta(:,end+1) = zeros(length(lamda),1);

% load R, f -------------------------------------------------
f = csvread('f.csv');
R = csvread('R2.csv');
    f(end) = [];
    R(end) = [];

theta = linspace(0,90*pi/180,10);
theta_interp = linspace(0,89*pi/180);
Lambertian = cos(theta_interp);
I_PD_interpolated = Lambertian;
% load Ipd(theta) angular profile --------------------------------
% prefix='angle d1';
all_files = dir;
file_names =cellstr('dummy');
for i=1:length(all_files)
    if(strfind(all_files(i).name,prefix)==1) 
        file_names = [file_names all_files(i).name];
    end
end
file_names(1) = []; % take out dummy
mat_names =cellstr('dummy');
    name_fix = strrep(file_names(1),' ','_');
    name_fix = strrep(name_fix,'.','_');
    mat_name = strcat('mat_',name_fix);
    mat_names = [mat_names mat_name];
    eval(char(strcat(mat_name,'=[];')));
    mat_names(1) = []; % take out dummy
   mat=[];
   fid = fopen(file_names{1});
    tline = fgets(fid);
    line_number = 1;
    while ischar(tline)
       line_number = line_number + 1;
       oline = strsplit(strtrim(tline));
       if ~isnan(str2double(oline(1))) % if the first entry in that line is a valid number, then u take that line and process and append to mat
           tmp_mat=zeros(1,length(oline));
           for j=1:length(oline)
               tmp_mat(j) = str2double(oline(j));
           end
           mat = cat(1,mat,tmp_mat);
       end
        tline = fgets(fid);
    end
    fclose(fid);
    
    I_PD_theta1 = mat(:,end);
    I_PD_theta1 = [I_PD_theta1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]; 
    I_PD_theta1(80:90) = I_PD_theta1(80:90)*1; %if edge is not darkened
    j = 1;
    I_PD_theta = [];
        
    for i = 1:10:91
        I_PD_theta (j) = mean(I_PD_theta1(i:i+9));
        j = j+1;
    end

    I_PD_interpolated = interp1(theta,I_PD_theta,theta_interp);
%     
 % plot Ipd(theta)-------------------------------------------
    subplot(223); hold all;
    plot(theta_interp*180/pi,I_PD_interpolated/max(I_PD_interpolated),'LineWidth',3.5);
    plot(theta_interp*180/pi,cos(theta_interp),'*');
    xlabel('theta'); ylabel('Ipd(A)');
    legend('Measured','Lambertian');
 %----------------------------------------------------------
    
 I_PD_interpolated = I_PD_interpolated/max(I_PD_interpolated);
 I_PD_interpolated';
 
   
for i=1:length(lamda)
    Intensity_interpolated(i,:) = interp1(theta,Intensity_lamda_theta(i,:),theta_interp);
end
A = zeros(1,length(theta_interp));
for i = 1:length(A)
    A(i) = (trapz(lamda,Intensity_interpolated(:,i)))^-1;
end

rho_lamda_theta = zeros(length(Intensity_interpolated), length(A));

 rho_lamda_theta = Intensity_interpolated .* repmat(A,length(lamda),1);

 
% CIE = csvread('CIE.csv');
% xbar = CIE(:,1);
% ybar = CIE(:,2);
% zbar = CIE(:,3);
% lamda_CIE = linspace(400,750,length(xbar))';
% 
% 
%    X = trapz(lamda_CIE,rho_lamda_theta(1:length(xbar),1).*xbar);
%    Y = trapz(lamda_CIE,rho_lamda_theta(1:length(xbar),1).*ybar);
%    Z = trapz(lamda_CIE,rho_lamda_theta(1:length(xbar),1).*zbar);
%    
%    X_CIE1931 = X/(X+Y+Z);
%    Y_CIE1931 = Y/(X+Y+Z);
%    CIE1931 = [X_CIE1931 Y_CIE1931]
%    d_CIE1931 = sqrt((X_CIE1931-.33)^2+(Y_CIE1931-.33)^2)
% str = sprintf('CIE 1931 X = %f , CIE 1931 Y = %f, CIE 1931 D = %f',X_CIE1931,Y_CIE1931,d_CIE1931);
% title(str);
 
figure; hold all;
peak_normalized_spectrum0 = rho_lamda_theta(:,1)/max(rho_lamda_theta(:,1));

for p = 1:9
    
    peak_normalized_spectrum0_ALL(:,p) = rho_lamda_thetaRAW(:,p)/max(rho_lamda_thetaRAW(:,p));
end

plot(lamda,peak_normalized_spectrum0,'LineWidth',3.5);
title('peak normalized forward spectrum');
hold off;
 
% plot spectrum -----------------------------------------
%  figure; hold all;
% i = linspace(1,100,10);
%  plot(lamda, rho_lamda_theta(:,i),'LineWidth',3);
%  legend('0','10','20','30','40','50','60','70','80');
% xlabel('lamda'); ylabel('normalized spectrum');
% hold off;
% figure; hold all;
% plot(lamda*1e9, rho_lamda_theta(:,1)/max(rho_lamda_theta(:,1)),'LineWidth',3);
% xlabel('Wavelength (nm)'); ylabel('Forward Spectrum (a.u)');
% hold off;
 % ------------------------------------------------------
 
G = dir('*.dat');
% plotEQE = figure;
% plotJV = figure;
plotPE = figure('units','normalized','outerposition',[0 0 1 1]);
% plotCE = figure;
% plotL = figure;

%Color Scheme-------------------------------------------
C = {'b','r','g','y','m','k','c',[.2 .5 1], [.2 .6 .2],[1 .6 .3],[1 .3 .8],[.3 .1 .91],[.6 .99 .05],[.99 .13 .1],[.1 .5 .3]};
if length(G)>length(C)
    for i = 1:length(G)-length(C)
    C = [C {rand(3,1)'}];
    end
end

%-------------------------------------------------------
for ii = 1:length(G)
    
    
     mat=[];
   fid = fopen(G(ii).name);
    tline = fgets(fid);
    line_number = 1;
   
    while ischar(tline)
       line_number = line_number + 1;
       oline = strsplit(strtrim(tline));
       if ~isnan(str2double(oline(1))) % if the first entry in that line is a valid number, then u take that line and process and append to mat
           tmp_mat=zeros(1,length(oline));
           for j=1:length(oline)
               tmp_mat(j) = str2double(oline(j));
           end
           mat = cat(1,mat,tmp_mat);
       end
        tline = fgets(fid);
    end
        fclose(fid);
 G(ii).name
        V = mat(:,1);
        J = mat(:,2)* 10; % mA/cm^2 to A/m^2
        I = J*Aoled;
        I_PD = mat(:,3);
    % JV--------------------------------------------------------------
         hold all;
        figure(plotPE);subplot(223);
        semilogy(V,J* (100^-2) * (1e3),'color',C{ii},'LineWidth',3);

        ylabel('Current Density (mA/cm^2)');
        xlabel('V');
        JVlegend{ii} = [G(ii).name]; % or whatever is appropriate

    % Ipd(theta) applied to JV-----------------------------------------    
        I_PD_V_theta = repmat(I_PD,1,length(I_PD_interpolated)).* repmat(I_PD_interpolated,length(V),1);

        
    % EQE-------------------------------------------------------------
        num2 = [];
        denom1 = [];
        EQE_num = [];
        num3 = [];
        num1 = [];

        num1 = I_PD_V_theta * (r^2) .* sin( repmat(theta_interp,length(V),1) );
        for i=1:length(theta_interp)
            num2(i) = trapz(lamda, rho_lamda_theta(:,i) .* lamda);
            denom1(i) = Apd* trapz(lamda, rho_lamda_theta(:,i) .* R );
        end

        for i = 1:length(V)
            EQE_num(i) = 2*pi* trapz(theta_interp, num1(i,:) .* num2 ./ denom1 );
        end

        EQE = 100*EQE_num' ./ ( (1240e-9)* I);
        EQE_ALL(:,ii) = EQE;
        J = J/10; %A/m^2 to mA/cm^2
        J_ALL(:,ii) = J;
        hold all;
        figure(plotSpec);subplot(224);
        semilogx(J,EQE,'color',C{ii},'LineWidth',3);
        ylim([0,50]);
        xlabel('current density (mA/cm^2)'); ylabel('EQE');
%         legendInfo{ii} = [G(ii).name]; % or whatever is appropriate

%     %PE----------------------------------------------------------------
        Radiance_W = zeros(1,length(I_PD));
        
        
        Radiance_lm = zeros(1,length(I_PD));
        Radiance_lm_0 = zeros(1,length(I_PD));
        num3 = [];
        PE_num = [];
        PE = [];
        CE = [];
        for i = 1:length(I_PD)
            Radiance_W (i) = (I_PD(i) * r^2)./(Apd * trapz(lamda,rho_lamda_theta(:,1).*R) ) ./ (Aoled*cos(0));
        end


        Radiance_lm = zeros(1,length(I_PD));
        for i=1:length(Radiance_W)
            Radiance_lm(i) = 683*Radiance_W(i)*trapz(lamda,rho_lamda_theta(:,1).*f);
        end
% Luminance----------------------------------------------------------
        Radiance_lm_0 = 683 * Radiance_W .* trapz(lamda,rho_lamda_theta(:,1).*f);
        Radiance_lm_0_ALL(:,ii) = Radiance_lm_0;
        hold all; figure(plotPE); subplot(224); 
        semilogy(V,Radiance_lm_0,'color',C{ii},'LineWidth',3.5); 
        xlabel('V'); ylabel('Forward Luminance (cd/m^2)')
% -------------------------------------------------------------------
        for i = 1:length(theta_interp)
            num3(i) = trapz(lamda, rho_lamda_theta(:,i) .* f );
        end

        for i = 1:length(V)
            PE_num(i) = trapz(theta_interp, num1(i,:) .* num3 ./ denom1 );
        end

        PE = 2*pi*683 * PE_num' ./ (I.*V);
        PE_ALL(:,ii) = PE;
        hold all;
        figure(plotPE); subplot(221);
        semilogx(Radiance_lm_0',PE,'color',C{ii},'LineWidth',3.5);
        xlabel('Forward Luminance (cd m^-2)');
        ylabel('Power Efficiency (lm/W)');
%         xlim([1 max(Radiance_lm_0)]);
      
     % CE ------------------------------------------------------------------
        
        CE = Radiance_lm_0 * Aoled ./ (I)';
          CE_ALL(:,ii) = CE;
        hold all;
        figure(plotPE); 
        subplot(222);
        semilogx(Radiance_lm_0',CE,'color',C{ii},'LineWidth',3.5);
        xlabel('Forward Luminance (cd m^-2'); ylabel('Current Efficacy (cd/A)');        
end
legend(JVlegend)

% PD Power---------------------------------------------------
for i = 1:length(theta_interp)
    Ppd_V_theta(:,i) = I_PD_V_theta(:,i) ./ trapz(lamda, rho_lamda_theta(:,i) .* R );
end
% figure; hold all;
% plot(theta_interp*180/pi,Ppd_V_theta(end,:));
% title('Total PD Power at 6.5V');
% xlabel('theta');



