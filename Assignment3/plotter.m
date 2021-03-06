%A = [1, 1, 1, 1, 1; 1, 1, 2, 3, 4; 1, 2, 1, 1, 1; 2, 1, 1, 2, 3];


%fileID = fopen('68_8.csv');
%A = textscan(fileID, '%d %d %f %f %f %f %f %f','delimiter', ',', 'EmptyValue', 0);
%A = dlmread('');
%   Reads this file, so make sure it is in the working folder!
A = csvread('64_8_2.csv');

%fclose(fileID);

n=size(A,2);
Psizes = [64,512,1024,2048,4096,8192];
Numproc = [8,16,32,64];
Sample = zeros(1,n-3);
Means = zeros(length(Numproc),n-3);
Variances = zeros(length(Numproc),n-3);
colours=['m','g','b','r','c','y'];

for P = Psizes
k=1;
for i = 1:size(A,1)
    if(A(i,1)==P)
        for j = Numproc
            if(A(i,2))== j
                Sample(k,:)=A(i,4:n);
                k=k+1;
            end
        end
    end
end

k=1;
i=1;
for j = Numproc
    Means(i,:) = mean(Sample(k:k+j-1,:),1);
    Variances(i,:) = var(Sample(k:k+j-1,:),1);
    k=k+j;
    i=i+1;
end
Means(:,1) = Means(:,1).*(Numproc.');

figure; hold on;
title(sprintf('Strong scaling for problem size %d',P))
xlabel('number of processors')
ylabel('time')
plot( Numproc, Means(:,1),'y')
for i = 2:size(Means,2)
    % shaded variance corridor uses external function 
    % https://de.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m?s_tid=srchtitle
    %plot(Numproc,Means(:,i))
    boundedline( Numproc, Means(:,i), sqrt(Variances(:,i)), colours(i-1) )
    alpha(.5)
end
grid on
%legend(h([1 3 5 7 9]),{'i/o','setup','compute','mpi','total'});
legend('i/o','setup var','setup','compute var','compute','mpi var','mpi','total var', 'total')
%legend('i/o','setup','compute','mpi','total')
% uses http://de.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline
vline(8), vline(16), vline(32), vline(64)

hold off;

end
