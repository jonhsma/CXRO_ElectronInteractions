function allocatedArray = arrayAllocator(N,varargin)
    % Check if there is a prototype
    if nargin >1 && isstruct(varargin{1})
        % Use the prototype
        allocatedArray(N) = varargin{1};
    else
        % Default declaration, which creates a prototype and throw it back
        % to the first case of this conditional
        % Post propagation values
        prototype.xyz_init  =   zeros([1 3]);
        prototype.xyz_final =   zeros([1 3]);
        prototype.xyz_delta =   zeros([1 3]);
        prototype.pathlen   =   -1.5;
        prototype.Eloss     =   -1.5;
        prototype.imfp      =   -1;
        prototype.deltaR    =   -1;   
        %%% parameters in scattering frame or San Francisco
        prototype.theta_SF      =   NaN;
        prototype.phi_SF        =   NaN;
        prototype.scattType     =   "Initialization";
        prototype.act           =   "Initialization";
        %%% parameters in lab frame
        prototype.theta_init  =   NaN;
        prototype.theta_final =   NaN;
        prototype.phi_init    =   NaN;
        prototype.phi_final   =   NaN;
        
        allocatedArray = arrayAllocator(N,prototype);
    end
end