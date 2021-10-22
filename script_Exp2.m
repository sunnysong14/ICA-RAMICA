% USAGE 
%   Script for the second experiment, i.e. bline image separation
% 
% Liyan for AAAI17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P_noise = [0:0.01:0.1, 0.15, 0.2];
StrICA_lst = {'infomax','FastICA','jade','ours'};
CRImg_lst = {'c','r'};
nS = 4;

% Caltech256 is 256x256
nRow = 64; %dimRVec
nCol = 256; %nSmp

% Special for RAMICA: if both c-&r- wise exit -> add one $StrICA_lst$
if sum(strcmpi(StrICA_lst,'ours')) && length(CRImg_lst)==2
    StrICA_lst{end+1} = 'ours'; %Add extra RAMICA at the end
end
nICA = length(StrICA_lst);  % rearrange.

p_noise = 0;
dir_data = ['..',filesep,'data',filesep,'Caltech256',filesep,'nS', num2str(nS), filesep];

seed_list = 100;

Amari = -ones(nICA, length(seed_list)); %init
for ic = 1 : nICA
    strICA = lower(StrICA_lst{ic});
    crImg = 'Not-used-for-this-ICA'; % overwritten by RAMICA.

    % special for RAMICA
    if strcmpi(strICA, 'ours')
        if ic==nICA && length(CRImg_lst)==2 %Second RAMICA
            crImg = CRImg_lst{2};
        else
            crImg = CRImg_lst{1};
        end
    end
   
    for ss = 1 : length(seed_list)
        seed = seed_list(ss);

        % load the source data in vector format
        flnmData = [dir_data, num2str(seed)];
        imgdata = load(flnmData); 
        scenstr = 'S_scen1';
        S_org = imgdata.(scenstr);

        % convert to data tensor: X_pn and TX_pn
        [X_pn,X,S,A,Noise_cov] = obtainICA_IS(S_org,p_noise,seed);
        %S: figure,for k=1:4,subplot(2,2,k),imshow(reshape(S(k, :),256,256),[]); end
        %X_pn: figure,for k=1:4,subplot(2,2,k),imshow(reshape(X_pn(k, :),256,256),[]); end

        % TX_pn: column-wise by default
        TX_pn = zeros(nS, nRow, nCol);
        for t = 1 : nCol
            TX_pn(:,:,t) = X_pn(:, ((t-1)*nRow+1):(t*nRow));
        end
        %Rest for row-wise
        if strcmpi(strICA, 'ours') && strcmpi(crImg, 'r')
            TX_pn_r = zeros(nS,nCol,nRow);
            for kk = 1 : nS
                TX_pn_r(kk,:,:)=squeeze(TX_pn(kk,:,:))'; %transpose
            end
            TX_pn = TX_pn_r; %replace
        end
        %figure,for k=1:4,subplot(2,2,k),imshow(squeeze(TX_pn(k,:,:)),[]);end

        % ICA method for source image recovery
        switch lower(strICA)                            
            case 'fastica' % 'pow3' (default)
                [RTS_est, ~, Ainv_est] = fastica(X);
            case 'jade'
                Ainv_est = jadeR(X);
                RTS_est = Ainv_est * X;
            case 'infomax'
                Ainv_est = infomaxICA(X);
                RTS_est = Ainv_est * X;
            case lower('ours')
                [Ainv_est, TS_est] = mJADE(TX_pn);
                RTS_est = TS_est;
            otherwise
                error('Error: undefined ICA method.');
        end

        % compute Amari error
        if size(Ainv_est,1) < size(A,1)
            amari = -1; 
        else
            P = Ainv_est * A;
            amari = Amari_index_ISA(P, ones(1,size(S,1)), 'uniform', 2);
        end
        
        Amari(ic, ss) = amari;
    end
    
end%ic

% ave PF across seeds
Ave_amari = -ones(size(Amari,1), 1);
Std_amari = -ones(size(Amari,1), 1);
for ic = 1 : size(Amari, 1)
    amari_arr = Amari(ic, :);
    pos_ = (amari_arr ~= -1);

    Ave_amari(ic) = mean(amari_arr(pos_));
    Std_amari(ic) = std(amari_arr(pos_));
end

% report the results
Ave_amari * 100


