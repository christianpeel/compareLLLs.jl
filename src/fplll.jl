
function initfplll()
#    const path_to_lib = "../fplll/"
    path_to_lib = "../fplll/"
    addHeaderDir(path_to_lib, kind=C_System)
    addHeaderDir("/usr/local/include", kind=C_System)
    Libdl.dlopen(path_to_lib * "/lib/libfplll.dylib", Libdl.RTLD_GLOBAL)

    cxxinclude("fplll.h")
    cxxinclude("fplll/main.h")
    cxx"""#include <iostream>"""
end

function showGMPinfoFromCxx()
    println("GMP info from Cxx:")
    print("  gmp_version is                ")
    icxx"""cout << gmp_version << endl;""";
    
    print("  GMP_LIMB_BITS is              ")
    icxx"""cout << GMP_LIMB_BITS << endl;""";

    println("  There is likely a way to check the path of the GMP library \n"*
            "    from within Cxx. For now we will assume that it is the \n"*
            "    same as in Julia.")
end

function showGMPinfoFromJulia()
    println("GMP info from Julia (GMP is used for BigInts):")
    println("  Base.GMP.VERSION is           $(Base.GMP.VERSION)  (compile-"*
          "time version)")
    if Base.GMP.VERSION==Base.GMP.version()
        println("  Base.GMP.version() is also    $(Base.GMP.version())  "*
              "(dynamic library version)")
    else
        warning("  This disagrees with Base.GMP.version(), the dynamic library "*
              "version, which is $(Base.GMP.version()).")
    end
    println("  Base.GMP.bits_per_limb() is   $(Base.GMP.bits_per_limb())")
    println("  Base.GMP.Limb is              $(Base.GMP.Limb) (datatype for Limb)")
    println("  path: $(dlpath(:libgmp))")
end

function showMPFRinfoFromCxx()
    println("MPFR info from Cxx:")
    print("  mpfr_version is           ")
    icxx"""cout << mpfr_version << endl;""";
    
    print("  MPFR_PREC_MIN is          ")
    icxx"""cout << MPFR_PREC_MIN << endl;""";
    
    print("  MPFR_PREC_MAX is          ")
    icxx"""cout << MPFR_PREC_MAX << endl;""";
    
    println("  There is likely a way to check the path of the GMP library \n"*
            "    from within Cxx. For now we will assume that it is the \n"*
            "    same as in Julia.")
end

function showMPFRinfoFromJulia()
    println("MPFR info from Julia (MPFR is used for BigFloats):")
    println("  Base.MPFR.VERSION is           $(Base.MPFR.VERSION)  (ignore-"*
          "this, it's the Julia Version)")
    println("  Base.MPFR.version() is         $(Base.MPFR.version())  "*
              "(dynamic library version)")
    println("  Base.MPFR.ROUNDING_MODE is     $(Base.MPFR.ROUNDING_MODE[])")
    println("  Base.MPFR.DEFAULT_PRECISION is $(Base.MPFR.DEFAULT_PRECISION[])")
    println("  precision(BigFloat) is         $(precision(BigFloat))")
    println("  path: $(dlpath(:libmpfr))")
end
