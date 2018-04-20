pub = rospublisher('/raw_vel');
msg = rosmessage(pub);

%% 1
sendVel(-0.1, 0.1, pub, msg);

t = timer;
t.StartDelay = 3;
t.TimerFcn = @(~,~) sendVel(0.0, 0.0, pub, msg);
start(t)


%% 2(run)
sendVel(0.2, 0.1, pub, msg);

t = timer;
t.StartDelay = 3;
t.TimerFcn = @(~,~) sendVel(0.0, 0.0, pub, msg);
start(t)

%% 2(process)
% experimental
load ex2.mat
dataset2.raw = dataset(:,[1 2 3]);

dataset2.velocity = diff(dataset2.raw);
dataset2.dt = dataset2.velocity(:,1);
dataset2.vl = dataset2.velocity(:,2);
dataset2.vr = dataset2.velocity(:,3);
dataset2.d = 0.245;
dataset2 = processNeato(dataset2);

% predicted
pred2.dt = dataset2.dt;
pred2.vl = 0.2 * ones(size(pred2.dt)) / 10;
pred2.vr = 0.1 * ones(size(pred2.dt)) / 10;
pred2.d = 0.245;
pred2 = processNeato(pred2);
plotNeato(dataset2, pred2);
legend("Experimental","Predicted")

%% 3
calc3.dt = 1;
calc3.t = 0:calc3.dt:10;
calc3.u = linspace(0,3.2,numel(calc3.t))';
calc3.r = [1/2*cos(calc3.u) 3/4*sin(calc3.u)];
calc3.v = diff(calc3.r)/calc3.dt;
calc3.vlin = sqrt(sum(calc3.v.^2,2));
calc3.T = [calc3.v./calc3.vlin zeros(size(calc3.vlin,1),1)];
calc3.w = cross(calc3.T(1:end-1,:), diff(calc3.T)/calc3.dt);
calc3.w = calc3.w(:,3);

%% 4 - bridge of death
% close all

input4.d = 0.25;
input4.dt = 0.1;
input4.t = (0:input4.dt:20)';
input4.u = linspace(0,2.9,numel(input4.t))';
input4.r = [0.3960*cos(2.65*(input4.u+1.4)), -0.99*sin(input4.u+1.4)];
input4.v = diff(input4.r)/input4.dt;
input4.vlin = sqrt(sum(input4.v.^2,2));
input4.T = [input4.v./input4.vlin zeros(size(input4.vlin,1),1)];

input4.w = cross(input4.T(1:end-1,:), diff(input4.T)/input4.dt);
input4.w = input4.w(:,3);

input4.vl = input4.vlin(1:end-1) - input4.w*input4.d/2;
input4.vr = input4.vlin(1:end-1) + input4.w*input4.d/2;

t = timer;
t.ExecutionMode = 'fixedDelay';
t.Period = input4.dt;
t.StartDelay = 5;
t.TasksToExecute = numel(input4.vl);
t.TimerFcn = @(myTimer,~) sendVel(input4.vl(myTimer.TasksExecuted), input4.vr(myTimer.TasksExecuted), pub, msg);
t.StopFcn = @(~,~) sendVel(0.0, 0.0, pub, msg);
% start(t)

figure
clf
hold on
plot(input4.r(:,1),input4.r(:,2))
quiver(input4.r(1:end-1,1),input4.r(1:end-1,2),input4.v(:,1),input4.v(:,2));
hold off
axis equal

%% 4 plot
input4.dt = input4.dt .* ones(size(input4.vl));

load bridge2.mat
dataset4.raw = dataset(:,[1 2 3]);
dataset4.velocity = diff(dataset4.raw);
dataset4.dt = dataset4.velocity(:,1);
dataset4.vl = dataset4.velocity(:,2);
dataset4.vr = dataset4.velocity(:,3);
dataset4.d = 0.245;

input4plot = input4;
input4plot .vl = input4plot.vl / 10;
input4plot .vr = input4plot.vr / 10;
input4plot = processNeato(input4plot );
plotNeato(input4plot ,processNeato(dataset4))
legend("Predicted", "Experimental")

%% functions
function sendVel(vl, vr, pub, msg)
    sprintf("VL: %.4d, VR: %.4d", vl, vr)
    msg.Data = [vl vr];
    send(pub, msg);
end

function out = processNeato(struct)
    % struct needs:
    % x
    % y
    % dt
    % d
    struct.w = (struct.vr - struct.vl)/struct.d;
    struct.v = mean([struct.vl struct.vr], 2);

    struct.theta = cumtrapz(struct.w);
    struct.drdt = struct.v.*[cos(struct.theta) sin(struct.theta)];

    struct.onValues = abs(struct.drdt(:,1))>0.001 | abs(struct.drdt(:,2))>0.001;
    struct.onRange = find(struct.onValues,1,'first'):find(struct.onValues,1,'last');
    struct.drdt = struct.drdt(struct.onRange,:);
    struct.dt = struct.dt(struct.onRange,:);

    struct.r = cumsum(struct.dt.*cumtrapz(struct.drdt));
    struct.x = struct.r(:,1);
    struct.y = struct.r(:,2);
    out = struct;
end

function plotNeato(varargin)
    figure
    axis equal
    hold on
    for i = 1:length(varargin)
        struct = varargin{i};
        plot(struct.x,struct.y,'*')
    end
    hold off
end