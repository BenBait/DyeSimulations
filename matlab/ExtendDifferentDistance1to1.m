%R6GtoRB ratio1:1
clear all;
close all;

L=43; % Size of the particle in nanometers
filename = ('../data/SiR6G-SiRB.xlsx');
p = xlsread(filename,'B2:B102'); %fluorescence intensity of SiRB
q = xlsread(filename,'C2:C102'); %fluorescence intensity of SiR6G
x = xlsread(filename,'A2:A102');
filename = ('../data/experimental_1to1.xlsx');
exp     = xlsread(filename, 'B2:B102');
dev_low = xlsread(filename, 'C2:B102');
dev_up  = xlsread(filename, 'D2:B102');

Fl_Table = table(x);
N_R6G=1267;% The Number of R6g molecules
N_RB=1525;% The Number of RB molecules

% creation of an aray of dyes
figure(1);

leg_str = {};
y = 2.0:0.1:3.0;
i_s = 1;
for LMin_ind=y
    leg_str{i_s} = strcat(num2str(LMin_ind),'nm');
    disp(leg_str(i_s));
    i_s = i_s + 1;
  
    clear B C
    disp(LMin_ind);

    % The size of matrix for placing the dye molecules
    % The maximum number of dye molecules is NxN in 2D and N^3 in 3D
    N=round(L/LMin_ind);

    for i = 1:N^3
        if i < N_R6G+1
            B(i)=1;  % R6G is +1
        else
            if i < N_R6G+N_RB+1
                B(i)=-1;  % RB is -1
            else
                B(i)=0;

            end
        end
    end

    num_trials = 20;
    
    ave = zeros(length(x), 1);
    for Fl_ind = 1:num_trials
        %Creation of matrix (a) with the random distribution of the dyes   
        NewIndex = randperm(numel(B));
        for i = 1:N^3
            C(i)=B(NewIndex(i));
        end
        
        a=reshape(C,N,N,N);
        % check for paired RB molecules
        % we don't have to check for R6G molecules because each R6G molecule
        % gets its own FRET PAIR SEARCH
        pd=zeros(N,N,N);

        % The search for neighbors and Identification of FRET
        % Calculation of FRET pairs for current min distance LMin_ind
        N_FRET = 0;

        buff = 0.5;
        b = LMin_ind - buff;
        c = LMin_ind + buff;
        % set the random number generator up for the normal dist
        rng(0, 'twister');

        R0=8.79; %fortser distance (nm) 
        QYR6G=0.95; %quantum yield of R6G
        QYRB=0.31; %quantum yield of RB
        AbR6G=0.005;%absorbance of R6G at 488nm (excitation wavelength)
        AbRB=0.001; %absorbance of RB at excitation wavelength 
        % FLUORENCE FOR THIS TRIAL; build it as we find FRET pairs
        Fl = zeros(length(x), 1);

        for i = 1:N
            for j = 1:N
                for z = 1:N
                    if a(i,j,z)== 1  % FRET is always assigned to R6G molecule

                       FRET=true; % True when no FRET pairs were associated with this particular R6G molecule
                       % rand returns on value from uniform distribution
                       % between 0 and 1
                       dist = (c-b)*rand + b;
                       % FRET efficiencey
                       E_FRET = (R0^6) / (R0^6 + dist^6);

                       % FRET PAIR SEARCH
                       % Do all these checks because each R6G molecule can be a part of only one FRET pair
                       % We also now flag in a separate array if a point is an RB molecule
                        if FRET && i-1 >0 && a(i-1,j,z)==-1 && ...
                           pd(i-1, j, z)==0 

                            Fl = Fl + AbR6G*QYRB*E_FRET*p;
                            N_FRET = N_FRET+1;
                            FRET=false;
                            pd(i-1, j, z) = 1;
                        end
                        if FRET && j-1 >0 && a(i,j-1,z)==-1 && ...
                           pd(i, j-1, z)==0

                            Fl = Fl + AbR6G*QYRB*E_FRET*p;
                            N_FRET=N_FRET+1;
                            FRET=false;
                            pd(i, j-1, z) = 1; 
                        end
                        if FRET && z-1>0 && a(i,j,z-1)==-1 && ...
                           pd(i, j, z-1)==0

                            Fl = Fl + AbR6G*QYRB*E_FRET*p;
                            N_FRET=N_FRET+1;
                            FRET=false;
                            pd(i, j, z-1) = 1; 
                        end 
                        if FRET && i+1 <= N && a(i+1,j,z)==-1 && ...
                           pd(i+1, j, z)==0

                            Fl = Fl + AbR6G*QYRB*E_FRET*p;
                            N_FRET=N_FRET+1;
                            FRET=false;
                            pd(i+1, j, z) = 1; 
                        end
                        if FRET && j+1 <= N  && a(i,j+1,z)==-1 && ...
                           pd(i, j+1, z)==0
                       
                            Fl = Fl + AbR6G*QYRB*E_FRET*p;
                            N_FRET=N_FRET+1;
                            FRET=false;
                            pd(i, j+1, z) = 1; 
                        end
                        if FRET && z+1 <=N && a(i,j,z+1)==-1 && ...
                           pd(i, j, z+1)==0
                       
                            Fl = Fl + AbR6G*QYRB*E_FRET*p;
                            N_FRET=N_FRET+1; 
                            FRET=false;
                            pd(i, j, z+1) = 1; 
                        end 
                       % check if diagonal is a neighbor
                   end 
                end 
            end
        end

        % We have been building this fluorescence during the trial as we find FRET pairs
        Fl = Fl + (N_R6G-N_FRET)*QYR6G*AbR6G*q + (N_RB-N_FRET)*QYRB*AbRB*p;

        Fl_Table.(strcat('F',int2str(LMin_ind*10) , 'l',int2str(Fl_ind))) = Fl;
        ave = ave + Fl;       
    end
    plot(x, ave/num_trials);
    hold on
end 

leg_str{i_s}   = "Experimental";
leg_str{i_s+1} = "Lower Deviation";
leg_str{i_s+2} = "Upper Deviation";

% Now get the x indices for the experimental data
x = xlsread(filename,'A2:A102');
plot(x, exp);
hold on
plot(x, dev_low);
hold on
plot(x, dev_up);

hold off
title("Intensity vs. Wavelength")
disp(leg_str)
legend(leg_str, 'Location', 'northwest')
writetable(Fl_Table,'../data/norm_dist_dat.csv')
