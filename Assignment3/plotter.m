%A = [1, 1, 1, 1, 1; 1, 1, 2, 3, 4; 1, 2, 1, 1, 1; 2, 1, 1, 2, 3];


%fileID = fopen('68_8.csv');
%A = textscan(fileID, '%d %d %f %f %f %f %f %f','delimiter', ',', 'EmptyValue', 0);
%A = dlmread('');
A = csvread('64_8.csv');

%fclose(fileID);

n=size(A,2);
Psizes = [64,512,1024,2048,4096,8192];
Numproc = [8,16,32,64];
Sample = zeros(1,n-3);
Means = zeros(length(Numproc),n-3);
Variances = zeros(length(Numproc),n-3);
colours=['y','m','g','b','r'];

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
for i = 1:size(Means,2)
    % shaded variance corridor uses external function 
    % https://de.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m?s_tid=srchtitle
    %plot(Numproc,Means(:,i))
    boundedline( Numproc, Means(:,i), sqrt(Variances(:,i)), ...
        colours(i) )
end
grid on
legend('i/o','-','setup','setup var','compute','compute var','mpi','mpi var','total', 'tatal var')
hold off;

end