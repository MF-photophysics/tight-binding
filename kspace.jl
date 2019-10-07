module Kspace

#require("mytools.jl")
include("mytools.jl")
using .Mytools

export kpts_path, kpts_grid, ewald, makennkp_rect
export modgrid


#in general, kpoints are represented as row vectors (covectors)
#  and real space points are represented as column vectors (contravariant).
#  So, if gvecs is a basis of kspace, gvecs[i,:] is the ith gvector.
#  If cell is a basis of rspace, cell[:,i] is the ith lattice vector.

#Note: gvecs = 2*pi*inv(cell)

#kpts[i,:] = ith kpt
function kpts_path(kpts,n)
  m = size(kpts,1)
  dim = size(kpts,2)
  ans = zeros(Float64,(n*(m-1),dim))
  p = kpts[1,:]
  k = 1
  for i=1:m-1
    dp = (kpts[i+1,:]-kpts[i,:])/n
    for j=1:n
      ans[k,:] = p
      p = p+dp
      k = k+1
    end
  end
  return ans
end

#for 2D
#gvecs[i,:] = ith gvector
function kpts_grid2d(gvecs,gridshape)
  ans = zeros(Float64,(prod(gridshape),2))
  l=1
  for i = 1:gridshape[1]
    for j = 1:gridshape[2]
      ans[l,:] = reshape([1.0*(i-1)/gridshape[1],1.0*(j-1)/gridshape[2]],(1,2)) * gvecs 
      l = l+1
    end
  end
  return ans
end

# for arbitrary dimensions
# size(gridshape,1) = size(gvecs,1)
function kpts_grid(gvecs,gridshape)
  anslen = prod(gridshape)
  ndims = size(gridshape,1)
  glen = size(gvecs,2)
  ans = zeros(Float64,(anslen,glen))
  icoords = zeros(Int,ndims)

  for ii = 0:anslen-1
    icoords = multiforhelp(gridshape,ii)
    tempvec = zeros(Float64,size(gvecs,2))
    for jj =1:ndims
      tempvec = tempvec + 1.0*icoords[jj]/gridshape[jj]*vec(gvecs[jj,:])
    end
    ans[ii+1,:] = tempvec
  end

  return ans
end

# rectangular grid neighbors (6 for 3d)
function makennkp_rect(gridshape)
  npts = prod(gridshape)
  ndims = size(gridshape,1)
  ans = zeros(Int,(npts*2*ndims,2+ndims))

  count = 1
  for ii=0:npts-1
    icoords = multiforhelp(gridshape,ii)
    for idim=1:ndims
      for j=-1:2:1
        i1coord = icoords + indicvec(ndims,idim,j)
        i2coord = [mod(i1coord[jj],gridshape[jj]) for jj=1:ndims]
        i2idx = flati(gridshape,i2coord)

        ans[count,1] = ii+1
        ans[count,2] = i2idx+1

        for jdim=1:ndims
          if(i1coord[jdim] <0)
            ans[count,2+jdim] = -1
          elseif(i1coord[jdim] >= gridshape[jdim])
            ans[count,2+jdim] = 1
          else
            ans[count,2+jdim] = 0
          end
        end

        count = count+1
      end
    end
  end

  return ans
end




#Ewald summation
#cell[:,i] is the ith lattice vector.
#qs = charges
#ts[:,i] = positions of charge i
#sig = erf parameter
#rcut, gcut = realspace/recipspace cutoffs
function ewald(cell,qs,ts,sig,rcut,gcut)
  shortrange = 0.0
  longrange = 0.0
  selfint = 0.0

  gvecs = 2*pi*inv(cell)
  vol = abs(det(cell))
  nq = size(qs,1)

  for ix=-rcut:rcut
    for iy=-rcut:rcut
      for iz=-rcut:rcut
        for i=1:nq
          for j=1:nq

            if(i==j && ix==0 && iy==0 && iz==0)
              continue
            end

            rdiff = norm(ts[:,i]-ts[:,j]+ cell*[ix,iy,iz])
            shortrange = shortrange + 0.5*qs[i]*qs[j]/rdiff * erfc(rdiff/sqrt(2.0)/sig)

          end
        end
      end
    end
  end

  for gx=-gcut:gcut
    for gy=-gcut:gcut
      for gz=-gcut:gcut

        if(gx==0 && gy==0 && gz==0)
          continue
        end

        gvec = vec(transpose([gx,gy,gz])*gvecs)

        strfac = 0.0

        for i=1:nq
          strfac = strfac+qs[i]*exp(1.0im*dot(gvec,ts[:,i]))
        end

        longrange = longrange + 2.0*pi/vol * exp(-0.5*sig^2 * norm(gvec)^2)/norm(gvec)^2 * norm(strfac)^2

      end
    end
  end

  for i=1:nq
    selfint = selfint + 1.0/sqrt(2.0*pi)/sig * qs[i]^2
  end

  return shortrange+longrange-selfint
end

#translates p2 by lattice vectors so that 
# it is as close as possible to p1
#cell[:,i] is the ith lattice vector.
function modgrid(cell,p1,p2)
  dpf = inv(cell)*(p2-p1)
  dpint = map(int,dpf)
  mf = multifor(fill(3,size(cell,1)))
  dist = norm(p2-p1)
  pans = p2
  for j=1:size(mf,1)
    p3 = p2 -cell*(dpint+vec(mf[j,:])-fill(1,size(cell,1)))
    dist1 =norm(p3-p1)
    if(dist1<dist)
      dist = dist1
      pans = p3
    end
  end

  return pans
end

end
