module Tb

#require("mytools.jl")
include("mytools.jl")
using .Mytools

export   tbhammaker3d

function tbhammaker3d(kpt, nhs, hopsdat; spars=false)
  #kpt in crystal coords (in units of recip latt. vecs
  #spas = constuct a sparse matrix

  if(spars)
    nhops = size(hopsdat,1)

    rcidxs = zeros(Int,nhops)
    vals = zeros(Complex64,nhops)
    rcidxs1 = Int[]
    vals1 = Complex128[]

    for ih in 1:nhops
      dcell = [round(Int64,hopsdat[ih,i]) for i = 1:3]
      from = round(Int64,hopsdat[ih,4])
      to = round(Int64,hopsdat[ih,5])
      val = hopsdat[ih,6] + 1.0im*hopsdat[ih,7]

      phase = exp(2.0im *pi* dot(dcell,kpt))
      #println(from,", ",to)

      rcidxs[ih] = flati1([nhs,nhs],[from-1,to-1])
      vals[ih] = phase*val
      if(dcell[1]!=0 || dcell[2]!=0 || dcell[3]!=0 || from!=to)

        push!(rcidxs1, flati1([nhs,nhs],[to-1,from-1]))
        push!(vals1,  conj(phase*val))
      end

    end

    rcidxs2 = cat(1,rcidxs,rcidxs1) #,dims=1)
    vals2 = cat(1,vals,vals1) #,dims=1)
    nhops2 = size(vals2,1)

    perm = sortperm(rcidxs2)
    rcidxss = rcidxs2[perm]
    valss = vals2[perm]

    rcidxu = zeros(Int,nhops2)
    rowidxu = zeros(Int,nhops2)
    colidxu = zeros(Int,nhops2)
    valsu = zeros(Complex128,nhops2)
    jj=1
    currnum = rcidxss[1]
    rcidxu[1] = currnum
    valsu[1] = valss[1]
    for ii=2:nhops2
      if(currnum<rcidxss[ii])
        currnum = rcidxss[ii]
        jj=jj+1
        rcidxu[jj] = currnum
      end
      valsu[jj] = valsu[jj] + valss[ii]
    end

    for ii=1:jj
      rc = multiforhelp1([nhs,nhs],rcidxu[ii])
      rowidxu[ii] = rc[1]+1
      colidxu[ii] = rc[2]+1
    end

    #println(sort(rowidxu[1:jj]))
    #println(colidxu[1:jj])
    #println(valsu[1:jj])

    ans = sparse(rowidxu[1:jj],colidxu[1:jj],valsu[1:jj])
  else

    ans = zeros(Complex128,(nhs,nhs))

    for ih in 1:size(hopsdat,1)
      dcell = [round(Int64,hopsdat[ih,i]) for i = 1:3]
      from = round(Int64,hopsdat[ih,4])
      to = round(Int64,hopsdat[ih,5])
      val = hopsdat[ih,6] + 1.0im*hopsdat[ih,7]

      phase = exp(2.0im *pi* dot(dcell,kpt))
      #println(from,", ",to)
      ans[from,to] = ans[from,to] + phase*val
      if(dcell[1]!=0 || dcell[2]!=0 || dcell[3]!=0 || from!=to)
        ans[to,from] = ans[to,from] + conj(phase*val)
      end

    end
  end
  return ans
end

end
