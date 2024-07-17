function read_hdf5_diagrams(fileCONN_fun,fileDISC_fun,fileCONN_as,fileDISC_as;Γ="g5",subtract_vev=false)
    conn_f  = h5read(fileCONN_fun,Γ)
    disc_f  = h5read(fileDISC_fun,Γ)
    conn_as = h5read(fileCONN_as ,Γ)
    disc_as = h5read(fileDISC_as ,Γ)
    # as a cross check get the number of stochastic sources from the hdf5 file
    T, L  = h5read(fileCONN_fun,"lattice")[1:2]
    
    #rescale disconnected pieces to match the common normalisation
    rescale_disc = (L^3)^2 /L^3
    # use the stochastic estimator of arXiv:1607.06654 eq. (14)
    disc_ff = unbiased_estimator(disc_f;rescale=rescale_disc,subtract_vev)
    disc_aa = unbiased_estimator(disc_as;rescale=rescale_disc,subtract_vev)
    disc_fa = unbiased_estimator(disc_f,disc_as;rescale=rescale_disc,subtract_vev)
    # rescale now the connected pieces appropriately
    rescale_connected!(conn_f ,L)
    rescale_connected!(conn_as,L)
    return conn_f, conn_as, disc_ff, disc_aa, disc_fa
end
# NOTE: Products of disconnected diagrams appear
#       For products of the same flavour I use the stochastic estimator of arXiv:1607.06654 eq. (14)
#       For products of different flavours we use a simple average over the sources. 
#       In both cases we need to take care of the relative time slices
function unbiased_estimator(discon1,discon2;rescale=1,subtract_vev=false)
    # (1) average over different hits
    # (2) average over all time separations
    # (3) normalize wrt. time and hit average
    nconf1, nhits1, T = size(discon1)
    nconf2, nhits2, T = size(discon2)
    nconf = min(nconf1,nconf2)
    timavg = zeros(eltype(discon1),(nconf,T))
    norm   = T*nhits1*nhits2
    if subtract_vev
        vev1 = vev_contribution(discon1)
        vev2 = vev_contribution(discon2)
        @inbounds for h in 1:nhits1, t in 1:T, conf in 1:nconf
            discon1[conf,h,t] = discon1[conf,h,t] - vev1[t]
        end
        @inbounds for h in 1:nhits1, t in 1:T, conf in 1:nconf
            discon2[conf,h,t] = discon2[conf,h,t] - vev2[t]
        end
    end
    for t in 1:T
        for t0 in 1:T
            Δt = mod(t-t0,T)
            @inbounds for hit1 in 1:nhits1, hit2 in 1:nhits2
                for conf in 1:nconf
                    loop1 = discon1[conf,hit1,t]
                    loop2 = discon2[conf,hit2,t0]
                    timavg[conf,Δt+1] += loop1*loop2
                end
            end
        end
    end
    @. timavg = rescale*timavg/norm
    # transpose the matrix so that it has the same layout as the connected pieces
    return timavg
end
function unbiased_estimator_threaded(discon1,discon2;rescale=1,subtract_vev=false)
    # (1) average over different hits
    # (2) average over all time separations
    # (3) normalize wrt. time and hit average
    nconf1, nhits1, T = size(discon1)
    nconf2, nhits2, T = size(discon2)
    if subtract_vev
        vev1 = vev_contribution(discon1)
        vev2 = vev_contribution(discon2)
    end
    
    # permute dimensions for better memory acces
    discon1 = permutedims(discon1,(2,3,1))
    discon2 = permutedims(discon2,(2,3,1))
    sum_loop2 = sum(discon2,dims=1)

    nconf = min(nconf1,nconf2)
    timavg = zeros(eltype(discon1),(T,nconf))
    norm   = T*nhits1*nhits2
    
    if subtract_vev
        @inbounds for conf in 1:nconf, h in 1:nhits1, t in 1:T
            discon1[h,t,conf] = discon1[h,t,conf] - vev1[t]
            discon2[h,t,conf] = discon2[h,t,conf] - vev2[t]
        end
    end
    @batch for conf in 1:nconf
        for t in 1:T
            for t0 in 1:T
                Δt = mod(t-t0,T)
                @inbounds for hit1 in 1:nhits1
                    loop1 = discon1[hit1,t,conf]
                    timavg[Δt+1,conf] += sum_loop2[1,t0,conf]*loop1
                end
            end
        end
    end
    @. timavg = rescale*timavg/norm
    return permutedims(timavg,(2,1))
end
function unbiased_estimator(discon;rescale=1,subtract_vev=false)
    # (1) average over different hits
    # (2) average over all time separations
    # (3) normalize wrt. time and hit average
    nconf, nhits, T = size(discon)
    timavg = zeros(eltype(discon),(nconf,T))
    norm   = T*div(nhits,2)^2
    hitsd2 = div(nhits,2)
    if subtract_vev
        vev = vev_contribution(discon)
        @inbounds for h in 1:nhits, t in 1:T, conf in 1:nconf
            discon[conf,h,t] = discon[conf,h,t] - vev[t]
        end
    end
    for t in 1:T
        for t0 in 1:T
            Δt = mod(t-t0,T)
            @inbounds for hit1 in 1:hitsd2, hit2 in hitsd2+1:nhits
                for conf in 1:nconf
                    loop1 = discon[conf,hit1,t]
                    loop2 = discon[conf,hit2,t0]
                    timavg[conf,Δt+1] += loop1*loop2
                end
            end
        end
    end
    @. timavg = rescale*timavg/norm
    return timavg
end
function vev_contribution(discon)
    vev = dropdims(mean(discon,dims=(1,2)),dims=(1,2)) 
    return vev
end
function rescale_connected!(corr,L)
    n1 = L^6/2 # from the norm used in HiReo
    n2 = L^3   # only keep a norm for the spatial volume
    @. corr  *= (n1/n2)
end
