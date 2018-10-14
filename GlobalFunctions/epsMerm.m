function eps=epsMerm(w,wp,q,qf,vf,gamma)
%%% usage: eps = epsMerm(w,wp,q,qf,vf,gamma)
%%% w: frequency in rad/s
%%% wp: center frequency in rad/s
%%% q: Momentum transfered in /m
%%% qf: Fermi Momentum in /m
%%% vf: Fermi velocity in m/s
%%% gamma: Oscillator gamma in /s [same unit as w]

tmpterm=epslind(w+sqrt(-1)*gamma,wp,q,qf,vf)-1;
tmpterm2=epslind(0,wp,q,qf,vf)-1;
eps=1+((1+sqrt(-1)*gamma./w).*tmpterm)./(1+(sqrt(-1)*gamma./w).*tmpterm./tmpterm2);

% % fprintf('epsMerm: q/qf = %.4f\n',q/qf);
% idx=find(isnan(eps)==1);
% if ~isempty(idx)
%     for i = 1:length(idx)
%         fprintf('epsMerm: eps=NaN at i = %d\n',i);
%     end
%     error('epsMerm: Resolve the eps=NaN issues!');
% end