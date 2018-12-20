%look into ButtonDown function
figure
w = waitforbuttonpress
if w == 0
    disp('Button click')
else
    disp('Key press')
end

%look into "CurrentPoint"
tryy = @fig ;

figure('ButtonDownFcn',tryy)

function fig(

% example
f=figure
x=1:10
y=15:24
plot(x,y)
get(f,'CurrentPoint')
a=gca
h=get(a,'CurrentPoint')
h(1,1:2) %this gives point of click in appropriate units

while ~strcmp(q1,'n')
    q1 = input('n?')
    display('not yet')
end
