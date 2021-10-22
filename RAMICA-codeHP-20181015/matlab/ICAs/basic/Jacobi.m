function Demix_mat = Jacobi(CM, Whiten_mat, nSource, nSensor, nSample, nCumMat, verbose)
% Usage:
%   Jacobi diagonalizing cumulant matrix. LY wrapped from the orginal JADE codes.
%
% INPUT ARGUMENTS
%       CM                  Cumulant matrices
%       MX_WHITEN           Whitening matrix
%       NSOURCE             #source
%       NSENSOR             #sensor
%          o                usually $nSource = nSensor$
%       NTIME               #(time step) or #(sampling)
%       NCUMMTX             #(Cumulant matrices)       
%       VERBOSE             Set to 0 for quiet operation
% 
% OUTPUT ARGUMENTS
%       Demix_mat           Demixing matrix such that $S_est=Demix_mat*X$
% 
% Liyan 12/23/2015 Wed. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Init
if 0  %<-- original set to 0
    % LY: using This seems to improve the performance of jade and jade_3rd.
    % 01-05-2016 Tuesday.
    %
    % Init by diagonalizing a *single* cumulant matrix.  
    % It seems to save some computation time `sometimes'.  
    % Not clear if initialization is really worth since Jacobi rotations are very efficient.  
    % On the other hand, it does not cost much...

	if verbose
        fprintf('jade -> Initialization of the diagonalization\n');
    end
	[V,D]	= eig(CM(:,1:nSource)); % Selectng a particular cumulant matrix.
	for u=1:nSource:nSource*nCumMat,         % Accordingly updating the cumulant set given the init
		CM(:,u:u+nSource-1) = CM(:,u:u+nSource-1)*V ; 
	end;
	CM	= V'*CM;
    
else	%% The dont-try-to-be-smart init
	V	= eye(nSource) ; % la rotation initiale
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%<-- Compute the square sum of diagonal entries of all cumulant matrices:
% Computing the initial value of the contrast
%------------------------------------------------------------------
Diag    = zeros(nSource,1) ;
On      = 0 ;
Range   = 1 : nSource ;
for im = 1 : nCumMat
    Diag  = diag(CM(:,Range));      %<-- diagonal entries of the im-th cumulant matrix
    On    = On + sum(Diag.*Diag);   %<-- sum of square of diagonal entries
    Range = Range + nSource ;
end
Off = sum(sum(CM.*CM)) - On ;       %<-- sum of square of off-diagonal entries

%<-- Liyan alternating codes for above:
%   o calculate $On$: sum of diagonal entries.
On2 = 0;
for i_ = 1 : nCumMat
    Block_ = CM(:, ((i_-1)*nSensor +1 : i_*nSensor));
    On2 = On2 + norm(diag(Block_))^2;
end
% 
Off2 = sum(sum(CM.*CM)) - On2;       %<-- sum of square of off-diagonal entries
%------------------------------------------------------------------

seuil	= 1.0e-6 / sqrt(nSample) ; % A statistically scaled threshold on `small' angles
encore	= 1;
sweep	= 0;                % sweep number
updates = 0;                % Total number of rotations
upds    = 0;                % Number of rotations in a given seep
g	= zeros(2,nCumMat);
gg	= zeros(2,2);
G	= zeros(2,2);
c	= 0 ;
s 	= 0 ;
ton	= 0 ;
toff	= 0 ;
theta	= 0 ;
Gain    = 0 ;

% Joint diagonalization proper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    fprintf('jade -> Contrast optimization by joint diagonalization\n'); 
end

while encore, encore = 0;   
    if verbose
        fprintf('jade -> Sweep #%3d',sweep); 
    end
    sweep = sweep + 1;
    upds  = 0 ; 
    Vkeep = V ;
    
    for p = 1 : nSource-1
        for q = p+1 : nSource
            Ip = p : nSource : nSource*nCumMat ;
            Iq = q : nSource : nSource*nCumMat ;

            %%% computation of Givens angle
            g	    = [ CM(p,Ip)-CM(q,Iq) ; CM(p,Iq)+CM(q,Ip) ];
            gg    = g*g';
            ton   = gg(1,1)-gg(2,2); 
            toff  = gg(1,2)+gg(2,1);
            theta = 0.5*atan2( toff , ton+sqrt(ton*ton+toff*toff) );
            Gain  = (sqrt(ton*ton+toff*toff) - ton) / 4 ;

            %%% Givens update
            if abs(theta) > seuil   % if Gain > 1.0e-3*On/nSource/nSource ,
                encore  = 1 ;
                upds    = upds    + 1;
                c	= cos(theta); 
                s	= sin(theta);
                G	= [ c -s ; s c ] ;  %<-- Givens rotation

                pair 		= [p; q] ;
                V(:,pair) 	= V(:, pair) * G;
                CM(pair,:)	= G' * CM(pair, :);
                CM(:,[Ip Iq]) 	= [ c*CM(:,Ip)+s*CM(:,Iq) -s*CM(:,Ip)+c*CM(:,Iq) ];
                %<-- CM is changing.

                On   = On  + Gain;
                Off  = Off - Gain;
                % fprintf('jade -> %3d %3d %12.8f\n',p,q,Off/On);

            end
        end
    end

    if verbose
      fprintf(' completed in %d rotations\n', upds); 
    end
    updates = updates + upds;
end

if verbose
    fprintf('jade -> Total of %d Givens rotations\n', updates); 
end

%%% A separating matrix
%   ===================
Demix_mat	= V' * Whiten_mat ;

%%% Permut the rows of the separating matrix Demix_mat to get the most energetic components first.
%%% Here the **signals** are normalized to unit variance.  Therefore, the sort is
%%% according to the norm of the columns of A = pinv(Demix_mat)

if verbose 
    fprintf('jade -> Sorting the components %s. \n', updates); 
end

A = pinv(Demix_mat);
[Ds,keys] = sort(sum(A.*A));
Demix_mat = Demix_mat(keys,:);
Demix_mat = Demix_mat(nSource:-1:1,:);      % Is this smart ?


% Signs are fixed by forcing the first column of Demix_mat to have non-negative entries.
if verbose
    fprintf('jade -> Fixing the signs %s.  \n', updates); 
end

b = Demix_mat(:,1);
signs = sign(sign(b)+0.1); % just a trick to deal with sign=0
Demix_mat = diag(signs) * Demix_mat;

end  %<-- END OF FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
