%edition five matrix process

%main program
% x is a vector 1*n
% mat is a matrix m*n,
% coherencematrix is a matrix,m*(nfft/2) which represent the coherence of x and each row in mat. 
%  pay attention to a problem,the dimension of coherencematrix is critical.   
% nfft should be the power of 2
function [coherencematrix]=matcoh1(mat1,mat2,fs,pfs,nw,k,sw)

    if nargin<7 sw = 5; end
    if nargin<6 k = 100; end
    if nargin<5 nw = 16; end

    samp = size(mat1,2);
    TR = 1/fs;
    nfft = 2^(ceil(log2(samp)));
    point = round(pfs*TR*nfft+1);
    point1 = max(point-floor(sw/2),1);
    point2 = min(point+floor(sw/2),samp);
    
    point = round(0.02*TR*nfft+1);
    point1 = max(point-floor(sw/2),1);
    
    points = [point1:point2];
    outv = dpss(size(mat1,2),nw,k);

    cft1_mat=single(geteigen_mat(mat1,k,TR,nfft,outv));
    if size(mat2,1)~=0
        cft2_mat=single(geteigen_mat(mat2,k,TR,nfft,outv));
    else
        cft2_mat=cft1_mat;
    end
    
    coherencematrix = 0;
    for i=1:length(points)
        coherencematrix = coherencematrix + abs(mtm_coh5(cft1_mat,cft2_mat,k,TR,nfft,outv,points(i)));
    end
    coherencematrix = coherencematrix/length(points);
    
end

%%%%%%%%%%%%%%%%%%%%5


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [msc]=mtm_coh5(cft1,cft2,k,TR,nfft,outv,p)

    nav=k;
    n1=nav+1;
    n2=nav+2;
    nfreqs=nfft/2+1;

    % kernel code
    if mode(nfreqs,2)==0
        blklof=1;
    else
        blklof=0;
    end
    trnrm=sqrt(2*nav-2);
    blkhif=nfreqs-1+blklof;

    %cancel cft calculation

    cft1 = squeeze(cft1(p,:,:));
    if size(cft1,1)==1
      cft1=cft1.';
    end
    cft2 = squeeze(cft2(p,:,:));
    if size(cft2,1)==1
      cft2=cft2.';
    end

    msc = abs(cft1.'*conj(cft2)).^2;
%     msc = msc./(sum(abs(cft1).^2)'*sum(abs(cft2).^2)+eps);
    msc = msc./(sum(abs(cft1).^2)'*sum(abs(cft2).^2)+0.001);
   
%     xsm2 = abs(cft1.'*conj(cft2)).^2;
%     jkmsc = sqrt(xsm2./(sum(abs(cft1).^2)'*sum(abs(cft2).^2)+eps));
%     NTmsc = trnrm*log((1.0+jkmsc)./(1.0-jkmsc+eps))/2;
%     msc = tanh(NTmsc/trnrm).^2;


    % xsm2= abs(sum(cft1(p,:).*conj(cft2(p,:)),2))^2;                     
    % jkmsc = xsm2/(sum(abs(cft1(p,:)).^2,2)*sum(abs(cft2(p,:)).^2,2));
    % NTmsc=trnrm*log((1.0+sqrt(jkmsc))/(1.0-sqrt(jkmsc)))/2;

       
end


function cft=geteigen_mat(x,k,TR,nfft,outv)
    outv=outv*sqrt(1/TR);
    nfreqs=nfft/2+1;
    %returnzero default is zero

    [sp,tp]=size(x);
    outv = repmat(outv(:),1,sp);
    x = repmat(x',k,1);
    taperdata=outv.*x;
    taperdata = reshape(taperdata,tp,k,sp);
    %add it
    pad_tap_dat=zeros(nfft,k,sp);
    pad_tap_dat(1:tp,1:k,1:sp)=taperdata(1:tp,1:k,1:sp);
    cft=fft(pad_tap_dat);
    %cft=fft(taperdata(1:nftt,:,:));
    cft=cft(1:nfreqs,:,:);
end
