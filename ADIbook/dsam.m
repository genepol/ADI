%dsam generates a 16x31 grid of diffusion coefficients
%on page 55 of my ADI book in blocks of 5x10.
D5 = ones(16,31);
    for ri = 1:5 %row ri
        for cj = 11:20 %col cj
            D5(ri,cj) = 4;
        end
        for cj = 21:31
            D5(ri,cj) = 36;
%            D5(ri,cj) = 1;
        end
    end
    for ri = 6:10 %row ri
        for cj = 1:10
            D5(ri,cj) = 16;
%            D5(ri,cj) = 4;
        end
        for cj = 11:20 %col cj
            D5(ri,cj) = 100;
%             D5(ri,cj) = 8;
        end
        for cj = 21:31
            D5(ri,cj) = 1600;
%            D5(ri,cj) = 4;
        end
    end
    for ri = 11:16 %row ri
        for cj = 1:10
            D5(ri,cj) = 9;
%            D5(ri,cj) = 1;
        end
        for cj = 11:20 %col cj
            D5(ri,cj) = 25;
%            D5(ri,cj) = 4;
        end
    end
nV = 16; nH = 31;
%D5 is the sample problem diffusion coefficient.
%The boundary condition is zero on columns = 0,nH+1 and rows = 0, nV+1.
sepd
return
