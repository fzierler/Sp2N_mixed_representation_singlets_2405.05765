function correlation_matrix(conn_f, conn_as, disc_ff, disc_aa, disc_fa;Nf_fun=2,Nf_as=3,disc_sign=+1)
    N1,T1 = size(conn_f)
    N2,T2 = size(conn_as)
    N3,T3 = size(disc_aa)
    N4,T4 = size(disc_ff)
    N5,T5 = size(disc_fa)
    N = minimum((N1,N2,N3,N4,N5))
    @assert T1 == T2 == T3 == T4 == T5 
    T = T1
    
    corr_matrix = zeros(eltype(conn_f),(2,2,N,T))

    @. corr_matrix[1,1,:,:] = conn_f[1:N,1:T]  - Nf_fun*disc_sign*disc_ff[1:N,1:T]
    @. corr_matrix[2,2,:,:] = conn_as[1:N,1:T] - Nf_as *disc_sign*disc_aa[1:N,1:T]
    @. corr_matrix[1,2,:,:] = sqrt(Nf_fun*Nf_as)*disc_sign*disc_fa[1:N,1:T]
    @. corr_matrix[2,1,:,:] = sqrt(Nf_fun*Nf_as)*disc_sign*disc_fa[1:N,1:T]
    return corr_matrix 
end
function _bin_correlator_matrix(corr;binsize=2)
    nop, nop, N, T = size(corr)
    corr_binned = zeros(eltype(corr), nop, nop, N÷binsize, T)
    for i in 1:N÷binsize
        for j in 1:binsize
            offset = (i-1)*binsize
            corr_binned[:,:,i,:] += corr[:,:,offset+j,:]/binsize
        end 
    end
    return corr_binned
end
function _swap_at_crossing(val,swap)
    N_eig, T = size(val)
    c = swap:T-swap+1
    @assert N_eig == 2
    val[1,c], val[2,c] = val[2,c], val[1,c] 
    return val
end
function eigenvalues(corr;swap=nothing,t0=1)
    eigvals_jk = eigenvalues_jackknife_samples(corr;t0)
    eigvals, Δeigvals = apply_jackknife(eigvals_jk;dims=2)
    if isnothing(swap)
        return eigvals, Δeigvals
    else
        _swap_at_crossing(eigvals, swap)
        _swap_at_crossing(Δeigvals,swap)
        return eigvals, Δeigvals
    end
end
function eigenvalues_eigenvectors(corr;swap=nothing,t0=1)
    eigvals_jk, eigvecs_jk = eigenvalues_eigenvectors_jackknife_samples(corr;t0)
    eigvals, Δeigvals = apply_jackknife(eigvals_jk;dims=2)
    eigvecs, Δeigvecs = apply_jackknife(eigvecs_jk;dims=3)
    if isnothing(swap)
        return eigvals, Δeigvals, eigvecs, Δeigvecs
    else
        _swap_at_crossing(eigvals, swap)
        _swap_at_crossing(Δeigvals,swap)
        # TODO: Sort eigenvectors accordingly
        # TODO: Make it work in this notation
        for i in 1:2
            eigvecs[i,:,:]  = _swap_at_crossing(eigvecs[i,:,:] ,swap)
            Δeigvecs[i,:,:] = _swap_at_crossing(Δeigvecs[i,:,:],swap)    
        end
        # Although: The matrix of eigenvectors should be time-independent
        return eigvals, Δeigvals, eigvecs, Δeigvecs
    end
end
function eigenvalues_jackknife_samples(corr;t0 = 1, imag_thresh = 1E-11)
    sample = delete1_resample(corr)
    nops, nconf, T = size(sample)[2:4]
    eigvals_jk = zeros(eltype(sample),(nops,nconf,T))
    for s in 1:nconf, t in 1:T
        # smaller values correspond to a faster decay, and thus correspond to a larger masses
        # use sortby to sort the eigenvalues by ascending eigen-energy of the meson state
        vals = eigen(sample[:,:,s,t],sample[:,:,s,t0]).values
        max_imag = maximum(imag.(vals)) 
        #max_imag > imag_thresh && @warn "imaginary part of $max_imag exceeds threshold of $imag_thresh"
        eigvals_jk[:,s,t] = real.(vals)
    end
    return eigvals_jk
end
function eigenvalues_eigenvectors_jackknife_samples(corr;t0 = 1)
    sample = delete1_resample(corr)
    nops, nconf, T = size(sample)[2:4]
    eigvals_jk = zeros(eltype(sample),(nops,nconf,T))
    eigvecs_jk = zeros(eltype(sample),(nops,nops,nconf,T))
    for s in 1:nconf, t in 1:T
        # smaller values correspond to a faster decay, and thus correspond to a larger masses
        # use sortby to sort the eigenvalues by ascending eigen-energy of the meson state
        sol = eigen(sample[:,:,s,t],sample[:,:,s,t0],sortby= x-> abs(x))
        eigvals_jk[:,s,t] = sol.values
        # I am unsure if the average over all eigenvectors is correct.
        for i in 1:2
            eigvecs_jk[:,i,s,t] = normalize(sol.vectors[:,i])
        end
    end
    return eigvals_jk, eigvecs_jk
end
# generate a resample of the original correlator matrix
function delete1_resample(corr_matrix)
    nops,nconf,T = size(corr_matrix)[2:end]
    samples = similar(corr_matrix)
    # temporary array for jackknife sampling
    tmp = zeros(eltype(corr_matrix),(nops,nops,nconf-1,T))
    for index in 1:nconf    
        for i in 1:index-1
            tmp[:,:,i,:] = corr_matrix[:,:,i,:]
        end
        for i in 1+index:nconf
            tmp[:,:,i-1,:] = corr_matrix[:,:,i,:]
        end
        # perform average after deleting one sample
        samples[:,:,index,:] = dropdims(mean(tmp,dims=3),dims=3)
    end
    return samples
end
# apply jackknife resampling along dimension dims
function apply_jackknife(obs::AbstractArray;dims::Integer)
    N  = size(obs)[dims]
    O  = dropdims(mean(obs;dims);dims)
    ΔO = dropdims(sqrt(N-1)*std(obs;dims,corrected=false);dims)
    return O, ΔO
end
function apply_jackknife(obs::AbstractVector)
    N  = length(obs)
    O  = mean(obs)
    ΔO = sqrt(N-1)*std(obs,corrected=false)
    return O, ΔO
end
# apply jackknife while ignoring NaNs
function nan_apply_jackknife(obs::AbstractArray;dims::Integer)
    N  = size(obs)[dims]
    O  = dropdims(nanmean(obs;dims);dims)
    ΔO = dropdims(sqrt(N-1)*nanstd(obs;dims,corrected=false);dims)
    return O, ΔO
end
function nan_apply_jackknife(obs::AbstractVector)
    N  = length(obs)
    O  = nanmean(obs)
    ΔO = sqrt(N-1)*nanstd(obs,corrected=false)
    return O, ΔO
end
function cov_jackknife_eigenvalues(evjk::AbstractArray)
    Nev, Nsamples, T = size(evjk) 
    covm = zeros(Nev,T,T)
    for N in 1:Nev
        c0 = (Nsamples-1)*cov(evjk[N,:,:],dims=1,corrected=false)
        covm[N,:,:] = Hermitian(c0)
    end
    return covm
end