%R6GtoRB ratio:1:0.1
clear all;
close all;
L=43; % Size of the particle in nanometers
filename = ('SiR6G-SiRB.xlsx');
p = xlsread(filename,'B2:B102'); %fluorescence intensity of SiRB
q = xlsread(filename,'C2:C102'); %fluorescence intensity of SiR6G
x = xlsread(filename,'A2:A102');

Fl_Table = table(x);
N_R6G=1274;% The Number of R6g molecules
N_RB=510;% The Number of RB molecules
% creation of an aray of dyes

for LMin_ind=2.0:0.1:3
    clear B C
    disp(LMin_ind);
    %LMin_ind=2.4;  % The size of a minimum cell (Minimum distance between dyes)in nanometers (Saquib used 3.5)
    N=round(L/LMin_ind); % The size of matrix for placing the dye molecules; The maximum number of dye molecules is NxN in 2D and N^3 in 3D
   
    for i = 1:N^3
        if i < N_R6G+1;
            B(i)=1;  % R6G is +1
        else
            if i < N_R6G+N_RB+1
                B(i)=-1;  % RB is -1
            else
                B(i)=0;

            end
        end
    end
    %disp(B)

    for Fl_ind = 1:20
        %Creation of matrix (a) with the random distribution of the dyes   
        NewIndex = randperm(numel(B));
        for i = 1:N^3
            C(i)=B(NewIndex(i));
        end
        
        a=reshape(C,N,N,N);

        % The search for neighbors and Identification of FRET
        % Calculation of FRET pairs N_FRET_LMin for min distance LMin
        N_FRET_LMin = 0;
        for i = 1:N
            for j = 1:N
                for z = 1:N
                 if a(i,j,z)== 1  % FRET is always assigned to R6G molecule, and only once
                    FRET=true; % True when no FRET pairs were associated with this particular R6G molecule
                    if FRET && i-1 >0 && a(i-1,j,z)==-1
                        N_FRET_LMin=N_FRET_LMin+1;
                        FRET=false;
                    end
                    if FRET && j-1 >0 && a(i,j-1,z)==-1
                        N_FRET_LMin=N_FRET_LMin+1;
                        FRET=false;
                    end
                    if FRET && z-1>0 && a(i,j,z-1)==-1
                        N_FRET_LMin=N_FRET_LMin+1;
                        FRET=false;
                    end 
                    if FRET && i+1 <= N && a(i+1,j,z)==-1
                        N_FRET_LMin=N_FRET_LMin+1;
                        FRET=false;
                    end
                    if FRET && j+1 <= N  && a(i,j+1,z)==-1
                        N_FRET_LMin=N_FRET_LMin+1;
                        FRET=false; 
                    end
                    if FRET && z+1 <=N && a(i,j,z+1)==-1
                       N_FRET_LMin=N_FRET_LMin+1; 
                       FRET=false;
                    end 
                 end 
                end 
            end
        end

        R0=8.79; %fortser distance (nm) 
        QYR6G=0.95; %quantum yield of R6G
        QYRB=0.31; %quantum yield of RB
        AbR6G=0.005;%absorbance of R6G at 488nm (excitation wavelength)
        AbRB=0.001; %absorbance of RB at excitation wavelength 
        Fl_Table.(strcat('F',int2str(LMin_ind*10) , 'l',int2str(Fl_ind))) = (N_R6G-N_FRET_LMin)*QYR6G*AbR6G*q+N_FRET_LMin*AbR6G*QYRB*R0^6/(R0^6+LMin_ind^6)*p+(N_RB-N_FRET_LMin)*QYRB*AbRB*p;
    end
end 

% figure
% plot(x, Fl)
writetable(Fl_Table,'dataTableFl0p1August12.csv')
%hold on
%xlswrite('data100.csv',[x,y])
%hold off
% figure
% plot(x, Fl_Table.(strcat('Fl',int2str(Fl_ind))))
    
