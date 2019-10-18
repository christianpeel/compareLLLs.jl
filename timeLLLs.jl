Pkg.add("Cxx")
Pkg.add("Nemo")
Pkg.add("LLLplus")
Pkg.add("Plots")

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
icxx"""o.method = LM_PROVED;""";    # see defs.h
# icxx"""o.int_type = ZT_LONG;""";      # use double  (Int64?). See defs.h line 197
# icxx"""o.float_type = FT_DOUBLE;""";  # use double  (Float64?). See defs.h line 206
function fplll()
    icxx"""status=lll_reduction(b,u,u_inv,o.delta,o.eta,o.method,o.float_type,o.precision,flags);"""
end


δ=.99;
ctx = Nemo.lll_ctx(δ, 0.51)
Nsamp = 6
Nvec = Int64.(round.(2 .^(1:.5:8)));
Nnum = length(Nvec);
times = zeros(3,Nnum,Nsamp);
for nx=1:Nnum
    N = Nvec[nx]
    
    icxx"""b.resize($N,$N);""";
    icxx"""u.resize($N,$N);""";

    S = Nemo.MatrixSpace(Nemo.ZZ, N,N);

    b = zeros(BigInt,N,N);
#    b = zeros(Int64,N,N);

    for sx = 1:Nsamp
        icxx"""b.gen_uniform(20);""";
        for ix = 1:N, jx = 1:N
            b[jx,ix] = icxx"""b[$(ix-1)][$(jx-1)].get_si();"""
        end
        bn = S(Nemo.fmpz.(b')); # fmpz is also GMP

        times[1,nx,sx] = @elapsed B,T = LLLplus.lll(b,δ)
        times[2,nx,sx] = @elapsed fplll()
        times[3,nx,sx] = @elapsed Bn, Tn = Nemo.lll_with_transform(bn,ctx);

    end
end
tave  = Statistics.mean(times,dims=3)[:,:,1]

plot(Nvec,tave',label=["LLLplus GMP","fplll GMP","Nemo GMP"],
                       linewidth=3,legend=:topleft)
xaxis!("log10(N)", :log10)
yaxis!("time", :log10)
savefig("timeVdim_20bits.png")

#logt = log10.(tave*1e6)./ repeat(log10.(Nvec'),2,1);
#plot(Nvec,logt')
#xaxis!("log10(N)", :log10)

# ratio = tave[2,:] ./ tave[3,:]
# plot(Nvec,ratio)
# xaxis!("log10(N)", :log10)
# yaxis!("ratio", :log10)




δ=.99;
ctx = Nemo.lll_ctx(δ, 0.51)
Nsamp = 8
N = 25
#Bvec = 4:4:60;
Bvec = 5:50:505;
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


