% c7p5.m 
% 
% symmetric recurrent nonlinear network

M=[0 -1 ; -1 0];  % connection matrix
h=[1 ; 1];        % input
iters=100;        % maximum number of iterations
eta=0.1;          % Euler step

clf;
subplot(1,2,1);
hold on;
beta=2.0;         % high gain for the activation function

                  % now step the starting conditions around the outside
                  % of a grid 
v1=0;
for v2=0:.25:3
  c7p1sub(M,beta,h,[v1 ; v2],eta,iters);
  drawnow;
end
v2=0;
for v1=0:.25:3
  c7p1sub(M,beta,h,[v1 ; v2],eta,iters);
  drawnow;
end
v1=3;
for v2=0:.25:3
  c7p1sub(M,beta,h,[v1 ; v2],eta,iters);
  drawnow;
end
v2=3;
for v1=0:.25:3
  c7p1sub(M,beta,h,[v1 ; v2],eta,iters);
  drawnow;
end

                  % now label the graph
axis('square'); 
grid;
title(sprintf('beta=%f',beta),'FontSize',20);
xlabel('v_1','FontSize',16); ylabel('v_2','FontSize',16);


subplot(1,2,2);
hold on;
beta=0.5;         % low gain for the activation function

                  % now step the starting conditions around the outside
                  % of a grid 
v1=0;
for v2=0:.25:3
  c7p1sub(M,beta,h,[v1 ; v2],eta,iters);
  drawnow;
end
v2=0;
for v1=0:.25:3
  c7p1sub(M,beta,h,[v1 ; v2],eta,iters);
  drawnow;
end
v1=3;
for v2=0:.25:3
  c7p1sub(M,beta,h,[v1 ; v2],eta,iters);
  drawnow;
end
v2=3;
for v1=0:.25:3
  c7p1sub(M,beta,h,[v1 ; v2],eta,iters);
  drawnow;
end

axis('square'); 
grid;
title(sprintf('beta=%f',beta),'FontSize',20);
xlabel('v_1','FontSize',16); ylabel('v_2','FontSize',16);

% my script________________________________________________________________

tau_exc = 10
tau_inh = 30

M=[(1.25-1)/tau_exc, (-1/tau_exc) ; 1/tau_inh, (-1-1)/tau_inh];  % connection matrix
h=[1 ; 1];        % input
iters=100;        % maximum number of iterations
eta=0.1;          % Euler step

clf;
subplot(1,2,1);
hold on;
beta=2.0;         % high gain for the activation function

                  % now step the starting conditions around the outside
                  % of a grid 

