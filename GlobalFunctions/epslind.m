function eps = epslind(w,wp,q,qf,vf)

e=1.6*1e-19;
hbar_J=6.626*1e-34/(2*pi);
% w=w-wp;
u=w./(q.*vf);
% if isnan(u)
%     error('epslind: u=NaN');
% end
z=q./(2.*qf);
% if isnan(z)
%     error('epslind: z=NaN');
% end
% % % 
% % % Chi2=e^2/(pi*hbar_J*vf);
% % % Chi2=Chi2./(8.85*1e-12); % to make Chi2 above unitless
% % % f1=1/2+1/(8.*z).*(g(z-u)+g(z+u));
% % % 
% % % for i = 1:length(u)
% % %     if z+u(i)<1
% % %         f2(i)=pi/2.*u(i);
% % %     else
% % %         if z+u(i)>1 & abs(z-u(i))<1
% % %             f2(i)=pi./(8.*z).*(1-(z-u(i)).^2);
% % %         else
% % %             if abs(z-u(i))>1
% % %                 f2(i)=0;
% % %             end
% % %         end
% % %     end
% % % end
% % % f2=0;
% % % % eps=1+Chi2./z.^2.*(f1+sqrt(-1).*f2);
% % % eps=1+1./z.^2.*f1;
% % % 
% % % idx=find(isnan(eps)==1);
% % % if ~isempty(idx)
% % %     for i = 1:length(idx)
% % %         fprintf('epslind: eps=NaN at i = %d\n',i);
% % %     end
% % %     error('epslind: Resolve the eps=NaN issues!');
% % % end

f=1/2+1./(8.*z).*(1-(z-u).^2).*log((z-u+1)./(z-u-1))+1./(8.*z)*(1-(z+u).^2).*log((z+u+1)./(z+u-1));

eps=1+3.*wp.^2./(q.^2.*vf.^2).*f;

if isnan(eps)
    dbg=1;
end

    function out=g(x)
        out=(1-x.^2).*log(abs((1+x)./(1-x)));
    end

end