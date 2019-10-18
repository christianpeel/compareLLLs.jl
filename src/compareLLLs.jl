"""

This is not a complete package. It's here only for discussion with others
who are interested in creating a Julia wrapper around fplll.

"""
module compareLLLs

using Cxx
using Libdl
import LLLplus
using Nemo

export
   timeVsSize,
   showGMPinfoFromCxx,
   showGMPinfoFromJulia,
   showMPFRinfoFromCxx,
   showMPFRinfoFromJulia


include("fplll.jl") # the "show" functions above, also initfplll
#initfplll()   # Initialize access to the fplll library.

@warn("compareLLLs.jl is not a Julia package, despite containing some of \n"*
      "         the files required for a package. At present, the most \n"*
      "         useful file is `timeLLLs.jl` in the main directory of the \n"*
      "         repo. Rather than typing \n"*
      "               julia> using compareLLLs \n"*
      "         we suggest direcly using the `timeLLLs.jl` script. ")

end # module compareLLLs
