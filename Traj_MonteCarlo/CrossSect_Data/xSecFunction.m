function y = xSecFunction(mode,angle,eInc,parameters)
    C_MODE_RESOLVED         =   parameters.C_MODE_RESOLVED;
    A_MODE_ANGLE_RESOLVED   =   parameters.A_MODE_ANGLE_RESOLVED;
    B_MODE_ANGLE_RESOLVED   =   parameters.B_MODE_ANGLE_RESOLVED;
    
    a = A_MODE_ANGLE_RESOLVED(angle,mode);
    b = B_MODE_ANGLE_RESOLVED(angle,mode);
    c = cc(b,mode,C_MODE_RESOLVED);
    
    y = a.*(eInc.^b).*exp(-c.*eInc)./(((b./c).^b).*exp(-b));
    
    function k = cc(b,m,param)
        k = param(1,m) + b.*param(2,m);
    end
    

end