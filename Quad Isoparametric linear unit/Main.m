close all;
clear all;
clc;


itype = [1,1,1,1];
ym = 1000000;
pr = 0.23;
thic = 1;
nint = 2;



Elementcoor = cell(1,4);
Elementcoor{1}=[1,0,0,1;1,1,0,0];
Elementcoor{2}=[0,-1,-1,0;1,1,0,0];
Elementcoor{3}=[0,-1,-1,0;0,0,-1,-1];
Elementcoor{4}=[1,0,0,1;0,0,-1,-1];
mark = zeros(4,8);
mark(1,:) = [1,2,3,4,7,8,5,6];
mark(2,:) = [3,4,0,0,0,0,7,8];
mark(3,:) = [7,8,0,0,0,0,11,12];
mark(4,:) = [5,6,7,8,11,12,9,10];

Elementstiff = cell(1,4);
for i = 1:length(mark(:,1))
    Elementstiff{i}=quads(i,itype(i),nint,thic,ym,pr,Elementcoor{i});
end

Structurestiff = zeros(12,12);
for i = 1:length(mark(:,1))
    for j = 1:length(mark(1,:))
        for k = 1:length(mark(1,:))
            if(mark(i,j)~=0 && mark(i,k)~=0)                
                Structurestiff(mark(i,j),mark(i,k)) = Elementstiff{i}(j,k);               
            else
            end
        
        end
        
    end
end
F = [0,-100000,0,0,0,0,0,-10000,0,0,0,0];
u = Structurestiff^-1*F';










figure(1);
for i = 1:length(Elementcoor)
    for j = 1:length(Elementcoor{1}(1,:))-1
        line(Elementcoor{i}(1,[j,j+1]),Elementcoor{i}(2,[j,j+1]),'Color','b','LineWidth',4);
        hold on
    end
    
end

axis([-2,2,-2,2]);

U  = cell(1,4);
for i = 1:length(mark(:,1))
    for j = 1:length(mark(1,:))
        if(mark(i,j)~=0)    
            U{i}(j) = u(mark(i,j));                 
        else
            U{i}(j) = 0; 
        end
        
    end
end

for i = 1:length(Elementcoor{i})
    temp = [];
    for j = 1:2:length(U{1}(:))-1
        temp(:,end+1) = U{i}([j,j+1])';
        
        
    end
    Elementcoor{i} = Elementcoor{i} + temp;

end


for i = 1:length(Elementcoor)
    for j = 1:length(Elementcoor{1}(1,:))-1
        line(Elementcoor{i}(1,[j,j+1]),Elementcoor{i}(2,[j,j+1]),'Color','r','LineWidth',1);
        hold on
    end
    
end

axis([-2,2,-2,2]);











