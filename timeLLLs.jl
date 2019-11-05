Pkg.add("Cxx")
Pkg.add("Nemo")
# Pkg.add("LLLplus")
Pkg.add("Plots")


########################################
# Int64 bases: plot of time vs length of smallest basis as delta is varied
########################################

using Cxx
using Libdl
import Nemo
import LLLplus
import Statistics
using Plots
gr();

# following assumes fplll is in a sibling dir. It could also be in a system directory.
const path_to_lib = "../fplll/"
addHeaderDir(path_to_lib, kind=C_System)
addHeaderDir("/usr/local/include", kind=C_System)
Libdl.dlopen(path_to_lib * "/lib/libfplll.dylib", Libdl.RTLD_GLOBAL)

cxxinclude("fplll.h")
cxxinclude("fplll/main.h")
cxx"""#include <iostream>"""


cxx"""ZZ_mat<long> b;"""
cxx"""ZZ_mat<long> u;"""
cxx"""ZZ_mat<long> u_inv;"""
cxx"""Options o;"""
cxx"""int status, flags;"""

icxx"""flags = LLL_DEFAULT;""";
icxx"""o.action = ACTION_LLL;""";
icxx"""o.eta = LLL_DEF_ETA;""";     # 0.51; see defs.h
icxx"""o.method = LM_FAST;""";    # see defs.h
icxx"""o.int_type = ZT_LONG;""";      # use double  (Int64?). See defs.h line 197
icxx"""o.float_type = FT_DOUBLE;""";  # use double  (Float64?). See defs.h line 206
# ?does generating or passing out unimodular slow the lll down?
function fplll()
    icxx"""status=lll_reduction(b,u,u_inv,o.delta,o.eta,o.method,o.float_type,o.precision,flags);"""
end

# I typically copy-and-pasted the code from here to the plot multiple times
# when changing things like the bit depth or largest dimension.

δv=[collect(.49:.05:.99); .995]
Nδ = length(δv);
bits = 25;
Nsamp = 8
N = 100
k = Int64(floor(N/2))
times = zeros(2,Nδ,Nsamp);
small = zeros(2,Nδ,Nsamp);
icxx"""b.resize($N,$N);""";
icxx"""u.resize($N,$N);""";
b = zeros(Int64,N,N);
bfp = zeros(Int64,N,N);

for nx=1:Nδ
    δ = δv[nx]
    ctx = Nemo.lll_ctx(δ, 0.51)
    
    icxx"""o.delta = $(δ);""";

    for sx = 1:Nsamp
        
        icxx"""b.gen_qary($k,$bits);""";  
        for ix = 1:N, jx = 1:N
            b[jx,ix] = icxx"""b[$(ix-1)][$(jx-1)].get_si();"""
        end

        times[1,nx,sx] = @elapsed B,T = LLLplus.lll(b,δ)
        small[1,nx,sx] = B[:,1]'*B[:,1]
        times[2,nx,sx] = @elapsed fplll()
        for ix = 1:N, jx = 1:N
            bfp[jx,ix] = icxx"""b[$(ix-1)][$(jx-1)].get_si();"""
        end
        small[2,nx,sx] = bfp[:,1]'*bfp[:,1]

    end

tave  = Statistics.mean(times,dims=3)[:,:,1];
save  = Statistics.mean(small,dims=3)[:,:,1];
display([δv'; tave])
display([δv'; save])
end

plot(save',tave',label=["LLLplus.lll Int64",
                       "fplll.lll_w_transform Int64"],
                       linewidth=3,legend=:topleft)
xaxis!("norm of smallest vector", :log10)
yaxis!("time", :log10)
title!("fplll's gen_qary(floor(N/2),$bits) used")

savefig("timeVsmallest_$(bits)bitsInt64.png")



########################################
# Int64 bases: plot of time vs dimenstion
########################################

using Cxx
using Libdl
import Nemo
import LLLplus
import Statistics
using Plots
gr();

# following assumes fplll is in a sibling dir. It could also be in a system directory.
const path_to_lib = "../fplll/"
addHeaderDir(path_to_lib, kind=C_System)
addHeaderDir("/usr/local/include", kind=C_System)
Libdl.dlopen(path_to_lib * "/lib/libfplll.dylib", Libdl.RTLD_GLOBAL)

cxxinclude("fplll.h")
cxxinclude("fplll/main.h")
cxx"""#include <iostream>"""


cxx"""ZZ_mat<long> b;"""
cxx"""ZZ_mat<long> u;"""
cxx"""ZZ_mat<long> u_inv;"""
cxx"""Options o;"""
cxx"""int status, flags;"""

icxx"""flags = LLL_DEFAULT;""";
icxx"""o.action = ACTION_LLL;""";
icxx"""o.delta = LLL_DEF_DELTA;"""; # 0.99
icxx"""o.eta = LLL_DEF_ETA;""";     # 0.51; see defs.h
icxx"""o.method = LM_FAST;""";    # see defs.h
icxx"""o.int_type = ZT_LONG;""";      # use double  (Int64?). See defs.h line 197
icxx"""o.float_type = FT_DOUBLE;""";  # use double  (Float64?). See defs.h line 206
# ?does generating or passing out unimodular slow the lll down?
function fplll()
    icxx"""status=lll_reduction(b,u,u_inv,o.delta,o.eta,o.method,o.float_type,o.precision,flags);"""
end

# I typically copy-and-pasted the code from here to the plot multiple times
# when changing things like the bit depth or largest dimension.

δ=.99;
bits = 35;
ctx = Nemo.lll_ctx(δ, 0.51)
Nsamp = 8
Nvec = Int64.(round.(2 .^(2:7)))
Nnum = length(Nvec);
times = zeros(2,Nnum,Nsamp);
small = zeros(2,Nnum,Nsamp);
for nx=1:Nnum
    N = Nvec[nx]
    k = Int64(floor(N/2))
    
    icxx"""b.resize($N,$N);""";
    icxx"""u.resize($N,$N);""";

    b = zeros(Int64,N,N);
    bfp = zeros(Int64,N,N);

    for sx = 1:Nsamp
        icxx"""b.gen_qary($k,$bits);""";  
        for ix = 1:N, jx = 1:N
            b[jx,ix] = icxx"""b[$(ix-1)][$(jx-1)].get_si();"""
        end

        times[1,nx,sx] = @elapsed B,T = LLLplus.lll(b,δ)
        small[1,nx,sx] = B[:,1]'*B[:,1]
        times[2,nx,sx] = @elapsed fplll()
        for ix = 1:N, jx = 1:N
            bfp[jx,ix] = icxx"""b[$(ix-1)][$(jx-1)].get_si();"""
        end
        small[2,nx,sx] = bfp[:,1]'*bfp[:,1]

    end
end

tave  = Statistics.mean(times,dims=3)[:,:,1];
save  = Statistics.mean(small,dims=3)[:,:,1];
[Nvec'; tave]
[Nvec'; save]

plot(Nvec,tave',label=["LLLplus.lll Int64",
                       "fplll.lll_w_transform Int64"],
                       linewidth=3,legend=:topleft)
xaxis!("log10(N)", :log10)
yaxis!("time", :log10)
title!("fplll's gen_qary(floor(N/2),$bits) used")

savefig("timeVdim_$(bits)bitsInt64.png")



########################################
# GMP bases: plot of time vs dimenstion
########################################



using Cxx
using Libdl
import Nemo
import LLLplus
import Statistics
using Plots
gr();

# following assumes fplll is in a sibling dir. It could also be in a system directory.
const path_to_lib = "../fplll/"
addHeaderDir(path_to_lib, kind=C_System)
addHeaderDir("/usr/local/include", kind=C_System)
Libdl.dlopen(path_to_lib * "/lib/libfplll.dylib", Libdl.RTLD_GLOBAL)

cxxinclude("fplll.h")
cxxinclude("fplll/main.h")
cxx"""#include <iostream>"""

cxx"""ZZ_mat<mpz_t> b;"""
cxx"""ZZ_mat<mpz_t> u;"""
cxx"""ZZ_mat<mpz_t> u_inv;"""
# cxx"""ZZ_mat<long> b;"""
# cxx"""ZZ_mat<long> u;"""
# cxx"""ZZ_mat<long> u_inv;"""
cxx"""Options o;"""
cxx"""int status, flags;"""

icxx"""flags = LLL_DEFAULT;""";
icxx"""o.action = ACTION_LLL;""";
icxx"""o.delta = LLL_DEF_DELTA;"""; # 0.99
icxx"""o.eta = LLL_DEF_ETA;""";     # 0.51; see defs.h
#icxx"""o.method = LM_WRAPPER;""";    # see defs.h
#icxx"""o.method = LM_HEURISTIC;""";    # see defs.h
icxx"""o.method = LM_PROVED;""";    # see defs.h
icxx"""o.int_type = ZT_MPZ;""";      # use GMP
icxx"""o.float_type = FT_MPFR;""";
# icxx"""o.method = LM_FAST;""";    # see defs.h
# icxx"""o.int_type = ZT_LONG;""";      # use double  (Int64?). See defs.h line 197
# icxx"""o.float_type = FT_DOUBLE;""";  # use double  (Float64?). See defs.h line 206
# does generating or passing out unimodular slow the lll down?
function fplll()
    icxx"""status=lll_reduction(b,u,u_inv,o.delta,o.eta,o.method,o.float_type,o.precision,flags);"""
end


δ=.99;
bits = 35;
ctx = Nemo.lll_ctx(δ, 0.51)
Nsamp = 8
Nvec = Int64.(round.(2 .^(2:6)))
Nnum = length(Nvec);
times = zeros(3,Nnum,Nsamp);
small = zeros(3,Nnum,Nsamp);
icxx"""o.precision = 100;""";
setprecision(BigFloat,100)
for nx=1:Nnum
    N = Nvec[nx]
    k = Int64(floor(N/2))
    
    icxx"""b.resize($N,$N);""";
    icxx"""u.resize($N,$N);""";

    S = Nemo.MatrixSpace(Nemo.ZZ, N,N);

    # b = zeros(Int64,N,N);
    # bfp = zeros(Int64,N,N);
    b = zeros(BigInt,N,N);
    bfp = zeros(BigInt,N,N);

    for sx = 1:Nsamp
        icxx"""b.gen_qary($k,$bits);""";  
        for ix = 1:N, jx = 1:N
            b[jx,ix] = icxx"""b[$(ix-1)][$(jx-1)].get_si();"""
        end
        bn = S(Nemo.fmpz.(b')); # fmpz is GMP

        times[1,nx,sx] = @elapsed B,T = LLLplus.lll(b,δ)
        small[1,nx,sx] = B[:,1]'*B[:,1]
        times[2,nx,sx] = @elapsed fplll()
        for ix = 1:N, jx = 1:N
            bfp[jx,ix] = icxx"""b[$(ix-1)][$(jx-1)].get_si();"""
        end
        small[2,nx,sx] = bfp[:,1]'*bfp[:,1]
        times[3,nx,sx] = @elapsed Bn, Tn = Nemo.lll_with_transform(bn,ctx);
        small[3,nx,sx] = Float64((Bn[1,:]*Bn[1,:]')[1,1])

    end
end

tave  = Statistics.mean(times,dims=3)[:,:,1];
save  = Statistics.mean(small,dims=3)[:,:,1];
[Nvec'; tave]
[Nvec'; save]

# plot(Nvec,tave',label=["LLLplus.lll Int64",
#                        "fplll.lll_w_transform Int64","Nemo.lll GMP"],
#                        linewidth=3,legend=:topleft)
plot(Nvec,tave',label=["LLLplus.lll GMP",
                       "fplll.lll_w_transform GMP","Nemo.lll GMP"],
                       linewidth=3,legend=:topleft)
xaxis!("log10(N)", :log10)
yaxis!("time", :log10)
title!("fplll's gen_qary(floor(N/2),$bits) used")

savefig("timeVdim_$(bits)bitsGMP.png")
#savefig("timeVdim_$(bits)bits.png")

#logt = log10.(tave*1e6)./ repeat(log10.(Nvec'),2,1);
#plot(Nvec,logt')
#xaxis!("log10(N)", :log10)

    bfp = zeros(Int64,N,N);
    for ix = 1:N, jx = 1:N
        bfp[jx,ix] = icxx"""b[$(ix-1)][$(jx-1)].get_si();"""
    end

# ratio = tave[2,:] ./ tave[3,:]
# plot(Nvec,ratio)
# xaxis!("log10(N)", :log10)
# yaxis!("ratio", :log10)




δ=.99;
ctx = Nemo.lll_ctx(δ, 0.51)
Nsamp = 8
N = 20
Bvec = 4:12:120;
#Bvec = 5:50:505;
Bnum = length(Bvec);
icxx"""b.resize($N,$N);""";
icxx"""u.resize($N,$N);""";
S = Nemo.MatrixSpace(Nemo.ZZ, N,N);
times = ones(3,Bnum,Nsamp);

#b = zeros(BigInt,N,N);
#b = zeros(Int64,N,N);

for bx=1:Bnum
    Bt = Bvec[bx]

    for sx = 1:Nsamp
        icxx"""b.gen_uniform($Bt);""";
        # for ix = 1:N, jx = 1:N
        #     b[jx,ix] = icxx"""b[$(ix-1)][$(jx-1)].get_si();"""
        # end
        # bn = S(Nemo.fmpz.(b')); # fmpz is also GMP
        setprecision(BigFloat,Bt+1)
        
        b = rand(0:big(2)^Bt,N,N)
        bn = S(Nemo.fmpz.(b')); # fmpz is also GMP

        times[1,bx,sx] = @elapsed B,T = LLLplus.lll(b,δ)
        times[2,bx,sx] = @elapsed fplll()
        times[3,bx,sx] = @elapsed Bn, Tn = Nemo.lll_with_transform(bn,ctx);
    end
end
tave  = Statistics.mean(times,dims=3)[:,:,1]


plot(Bvec,tave',label=["LLLplus GMP","fplll GMP","Nemo GMP"],
                       linewidth=3,legend=:topleft)
xaxis!("B (Number of bits)")
yaxis!("time", :log10)
title!("$(N)x$(N) matrices were used")

savefig("timeVnbits_Neq25.png")


