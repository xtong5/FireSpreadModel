function [cmap]=buildcmap(colors)
% color map

if nargin<1
    colors='wrgbcmyk';
end
if ~ischar(colors)
    error(['Error! colors must be a variable of type char with '...
        'color-names, such as ''r'', ''g'', etc., ']);
end
ncolors=length(colors)-1;
bins=round(255/ncolors);
vec=zeros(300,3);
switch colors(1)
    case 'w'
        vec(1,:)=1;
    case 'r'
        vec(1,:)=[1 0 0];
    case 'g'
        vec(1,:)=[0 1 0];
    case 'b'
        vec(1,:)=[0 0 1];
    case 'c'
        vec(1,:)=[0 1 1];
    case 'm'
        vec(1,:)=[1 0 1];
    case 'y'
        vec(1,:)=[1 1 0];
    case 'k'
        vec(1,:)=[0 0 0];
end
for i=1:ncolors
 beG=(i-1)*bins+1;
 enD=i*bins+1; 
 switch colors(i+1)
     case 'w'
         vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';
     case 'r'
         vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';
     case 'g'
         vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';
     case 'b'         
         vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';
     case 'c'
         vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';
     case 'm'
         vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';
     case 'y'
         vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';
     case 'k'
         vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';
 end
end
cmap=vec(1:bins*ncolors,:);
end
