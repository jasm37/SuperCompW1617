
A = csvread('64_8_2.csv');

n=size(A,2);
Psizes = [64,512,1024,2048,4096,8192];
Numproc = [8,16,32,64];
Sample = zeros(1,n-3);
Means = zeros(length(Psizes),n-3);
Variances = zeros(length(Psizes),n-3);
colours=['m','g','b','r','c','y'];

for N = Numproc
k=1;
for i = 1:size(A,1)
    if(A(i,2)==N)
        for j = Psizes 
            if(A(i,1))== j
                Sample(k,:)=A(i,4:n);
                k=k+1;
            end
        end
    end
end

k=1;
i=1;
for j = Psizes 
    Means(i,:) = mean(Sample(k:k+N-1,:),1);
    Variances(i,:) = var(Sample(k:k+N-1,:),1);
    k=k+N;
    i=i+1;
end
Means(:,1) = Means(:,1).*N;

%------------------------------------------------------
figure; hold on;
title(sprintf('Strong scaling for problem size %d',N))
xlabel('number of processors')
ylabel('time')
plot(Psizes, Means(:,1),'y')
for i = 2:size(Means,2)
    % shaded variance corridor uses external function 
    % https://de.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m?s_tid=srchtitle
    %plot(Numproc,Means(:,i))
    boundedline( Psizes, Means(:,i), sqrt(Variances(:,i)), colours(i-1) )
    alpha(.5)
end
grid on
%legend(h([1 3 5 7 9]),{'i/o','setup','compute','mpi','total'});
legend('i/o','setup var','setup','compute var','compute','mpi var','mpi','total var', 'total')
% uses http://de.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline
vline(64), vline(512), vline(1024), vline(2048), vline(4096), vline(8192)
hold off;

end