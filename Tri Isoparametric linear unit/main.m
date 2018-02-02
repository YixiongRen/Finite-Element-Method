close all;
clear all;
clc;


itype = [1,1,1,1];
ym = 1000000;
pr = 0.23;
thic = 1;
nint = 3;



Elementcoor = cell(1,4);
Elementcoor{1}=[2,1,2;2,1,0];
Elementcoor{2}=[2,0,1;2,2,1];
Elementcoor{3}=[1,0,0;1,2,0];
Elementcoor{4}=[2,1,0;0,1,0];
mark = zeros(4,6);
mark(1,:) = [1,2,3,4,5,6];
mark(2,:) = [1,2,0,0,3,4];
mark(3,:) = [3,4,0,0,0,0];
mark(4,:) = [5,6,3,4,0,0];

Elementstiff = cell(1,4);
for i = 1:length(mark(:,1))
    Elementstiff{i}=tris(i,itype(i),nint,thic,ym,pr,Elementcoor{i});
end

Structurestiff = zeros(6,6);
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
F = [0,-100000,0,0,0,0];
u = Structurestiff^-1*F';










figure(1);
for i = 1:length(Elementcoor)
    for j = 1:length(Elementcoor{1}(1,:))-1
        line(Elementcoor{i}(1,[j,j+1]),Elementcoor{i}(2,[j,j+1]),'Color','b','LineWidth',4);
        hold on
    end
    
end

axis([-1,3,-1,3]);

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

axis([-1,3,-1,3]);