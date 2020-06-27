%init_guess
%generate an initial guess of the loading vector and factor matrices
%
%By Arie Yeredor, Nov. 2018, revised May 2020

minit=3;     %initialization mode
switch minit
    case 1  %true values
        if FF==F
            AA=Aall;
            ll=lam;
        else
            error('can''t initialize as true values when FF is not F');
        end
    case 2  %random values    
        rng(123)        %for repeatability; comment out for "true" randomness    
        AA=cell(1,N);
        for n=1:N
            Pn=rand(Ivec(n),FF);
            Pn=Pn./sum(Pn);
            AA{n}=Pn;
        end
        ll=rand(FF,1);
        ll=ll/sum(ll);
    case 3  %by the FF most probable according to Xemp
        AA=cell(1,N);
        for n=1:N
            AA{n}=0.001*ones(Ivec(n),FF);
        end
        ll=zeros(FF,1);
        intoX=cell(1,N);
        [sx ix]=sort(Xemp(:),'descend');
        for f=1:FF
            [intoX{:}]=ind2sub(Ivec,ix(f));
            for n=1:N
                AA{n}(intoX{n},f)=1;
                AA{n}=AA{n}./sum(AA{n});
            end
            ll(f)=sx(f);
        end
        ll=ll/sum(ll);
end
