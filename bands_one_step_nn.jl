module Bands


#require("htools.jl")
include("htools.jl")
using .Htools

#require("tb.jl")
include("tb.jl")
using .Tb

#require("mytools.jl")
include("mytools.jl")
using .Mytools

#require("kspace.jl")
include("kspace.jl")
using .Kspace

using Dates

export gap

function gap(it,nx,dat)

ncells = nx^3
nt = Int64(round(size(dat)[1]/(ncells*10)))
nat = ncells*12
count = 0


#epbs = -2.473
#epbp = 6.017
#eip = 2.518

#tsp = 1.0
#tpp = 1.8
#tpi = -0.4
socpb = 0.56
soci = 0.31

nnpb = 0.0838
nni = 0.155 #across a Pb
nni2 = 0.0450 #not across a Pb
use_pbpb_hop = true #hopping between adjacent Pb
use_ii_hop1 = true #hopping between I separated by a Pb
use_ii_hop2 = true # hopping between I not separated by a Pb
  println(it)

  hopsdat = [] #{}

 
  for ix=1:nx
    for iy=1:nx
      for iz=1:nx

        icell = iz+nx*(iy-1)+nx*nx*(ix-1)

        icellxm = iz+nx*(iy-1)+nx*nx*(mod1(ix-1,nx)-1)
        icellym = iz+nx*(mod1(iy-1,nx)-1)+nx*nx*(ix-1)
        icellzm = mod1(iz-1,nx)+nx*(iy-1)+nx*nx*(ix-1)
        icellms = [icellxm,icellym,icellzm]

        iblock = (it-1)*ncells+icell
        #println(iblock)

        epbs = dat[10*(iblock-1)+1,1]
        epbpx = dat[10*(iblock-1)+1,2]
        epbpy = dat[10*(iblock-1)+1,3]
        epbpz = dat[10*(iblock-1)+1,4]
        epbp = [epbpx,epbpy,epbpz]
        #println(epbp)
        eixpx = dat[10*(iblock-1)+2,2]
        eixpy = dat[10*(iblock-1)+2,3]
        eixpz = dat[10*(iblock-1)+2,4]
        eiypx = dat[10*(iblock-1)+3,2]
        eiypy = dat[10*(iblock-1)+3,3]
        eiypz = dat[10*(iblock-1)+3,4]
        eizpx = dat[10*(iblock-1)+4,2]
        eizpy = dat[10*(iblock-1)+4,3]
        eizpz = dat[10*(iblock-1)+4,4]
        eip = [eixpx, eixpy, eixpz, eiypx, eiypy, eiypz, eizpx, eizpy, eizpz ]

        #println(eip)

        tspsigxp = dat[10*(iblock-1)+5,1] 
        tspsigyp = dat[10*(iblock-1)+6,1] 
        tspsigzp = dat[10*(iblock-1)+7,1] 
        tspsigxm = dat[10*(iblock-1)+8,1] 
        tspsigym = dat[10*(iblock-1)+9,1] 
        tspsigzm = dat[10*(iblock-1)+10,1] 

        tppsigxp = dat[10*(iblock-1)+5,2] 
        tppsigyp = dat[10*(iblock-1)+6,2] 
        tppsigzp = dat[10*(iblock-1)+7,2] 
        tppsigxm = dat[10*(iblock-1)+8,2] 
        tppsigym = dat[10*(iblock-1)+9,2] 
        tppsigzm = dat[10*(iblock-1)+10,2] 

        tpppip = dat[10*(iblock-1)+5:10*(iblock-1)+7,3:4] 
        tpppim = dat[10*(iblock-1)+8:10*(iblock-1)+10,3:4] 


        #onsite
        for i=1:1
          push!(hopsdat,[0,0,0,(icell-1)*26+i,(icell-1)*26+i,epbs,0.0])
          push!(hopsdat,[0,0,0,(icell-1)*26+i+13,(icell-1)*26+i+13,epbs,0.0])
        end

        for i=2:4
          val = epbp[i-1]
          push!(hopsdat,[0,0,0,(icell-1)*26+i,(icell-1)*26+i,val,0.0])
          push!(hopsdat,[0,0,0,(icell-1)*26+i+13,(icell-1)*26+i+13,val,0.0])
        end

        for i=5:13
          val = eip[i-4]
          push!(hopsdat,[0,0,0,(icell-1)*26+i,(icell-1)*26+i,val,0.0])
          push!(hopsdat,[0,0,0,(icell-1)*26+i+13,(icell-1)*26+i+13,val,0.0])
        end
 
        #sp
        push!(hopsdat,[0,0,0,(icell-1)*26+1,(icell-1)*26+5,tspsigxp,0.0])
        push!(hopsdat,[0,0,0,(icell-1)*26+1,(icell-1)*26+9,tspsigyp,0.0])
        push!(hopsdat,[0,0,0,(icell-1)*26+1,(icell-1)*26+13,tspsigzp,0.0])

        push!(hopsdat,[0,0,0,(icell-1)*26+1+13,(icell-1)*26+5+13,tspsigxp,0.0])
        push!(hopsdat,[0,0,0,(icell-1)*26+1+13,(icell-1)*26+9+13,tspsigyp,0.0])
        push!(hopsdat,[0,0,0,(icell-1)*26+1+13,(icell-1)*26+13+13,tspsigzp,0.0])

        dx = (ix==1 ? -1 : 0)
        dy = (iy==1 ? -1 : 0)
        dz = (iz==1 ? -1 : 0)
        ds = [dx,dy,dz]

        push!(hopsdat,[dx,0,0,(icell-1)*26+1,(icellxm-1)*26+5,-tspsigxm,0.0])
        push!(hopsdat,[0,dy,0,(icell-1)*26+1,(icellym-1)*26+9,-tspsigym,0.0])
        push!(hopsdat,[0,0,dz,(icell-1)*26+1,(icellzm-1)*26+13,-tspsigzm,0.0])

        push!(hopsdat,[dx,0,0,(icell-1)*26+1+13,(icellxm-1)*26+5+13,-tspsigxm,0.0])
        push!(hopsdat,[0,dy,0,(icell-1)*26+1+13,(icellym-1)*26+9+13,-tspsigym,0.0])
        push!(hopsdat,[0,0,dz,(icell-1)*26+1+13,(icellzm-1)*26+13+13,-tspsigzm,0.0])

        #ppsi
        push!(hopsdat,[0,0,0,(icell-1)*26+2,(icell-1)*26+5,tppsigxp,0.0])
        push!(hopsdat,[0,0,0,(icell-1)*26+3,(icell-1)*26+9,tppsigyp,0.0])
        push!(hopsdat,[0,0,0,(icell-1)*26+4,(icell-1)*26+13,tppsigzp,0.0])

        push!(hopsdat,[0,0,0,(icell-1)*26+2+13,(icell-1)*26+5+13,tppsigxp,0.0])
        push!(hopsdat,[0,0,0,(icell-1)*26+3+13,(icell-1)*26+9+13,tppsigyp,0.0])
        push!(hopsdat,[0,0,0,(icell-1)*26+4+13,(icell-1)*26+13+13,tppsigzp,0.0])

        #println(bondxm)
        #println(-ppsi(norm(bondxm)))

        push!(hopsdat,[dx,0,0,(icell-1)*26+2,(icellxm-1)*26+5,tppsigxm,0.0])
        push!(hopsdat,[0,dy,0,(icell-1)*26+3,(icellym-1)*26+9,tppsigym,0.0])
        push!(hopsdat,[0,0,dz,(icell-1)*26+4,(icellzm-1)*26+13,tppsigzm,0.0])

        push!(hopsdat,[dx,0,0,(icell-1)*26+2+13,(icellxm-1)*26+5+13,tppsigxm,0.0])
        push!(hopsdat,[0,dy,0,(icell-1)*26+3+13,(icellym-1)*26+9+13,tppsigym,0.0])
        push!(hopsdat,[0,0,dz,(icell-1)*26+4+13,(icellzm-1)*26+13+13,tppsigzm,0.0])

        #pppi
        for i=0:2
          js = deleten([0,1,2],i+1)
          for j1=1:2
            j = js[j1]
            push!(hopsdat,cat(1,[0,0,0], [(icell-1)*26+2+j,(icell-1)*26+5+3*i+j,tpppip[i+1,j1],0.0]))
            push!(hopsdat,cat(1,indicvec(3,i+1,1)*ds[i+1], [(icell-1)*26+2+j,(icellms[i+1]-1)*26+5+3*i+j,tpppim[i+1,j1],0.0]))

            push!(hopsdat,cat(1,[0,0,0], [(icell-1)*26+2+j+13,(icell-1)*26+5+3*i+j+13,tpppip[i+1,j1],0.0]))
            push!(hopsdat,cat(1,indicvec(3,i+1,1)*ds[i+1], [(icell-1)*26+2+j+13,(icellms[i+1]-1)*26+5+3*i+j+13,tpppim[i+1,j1],0.0]))
          end
        end

        #soc

        push!(hopsdat,[0,0,0,(icell-1)*26+4,(icell-1)*26+2+13,-socpb,0.0])
        push!(hopsdat,[0,0,0,(icell-1)*26+4,(icell-1)*26+3+13,0.0,socpb])
        push!(hopsdat,[0,0,0,(icell-1)*26+2,(icell-1)*26+4+13,socpb,0.0])
        push!(hopsdat,[0,0,0,(icell-1)*26+3,(icell-1)*26+4+13,0.0,-socpb])
        push!(hopsdat,[0,0,0,(icell-1)*26+2,(icell-1)*26+3,0.0,-socpb])
        push!(hopsdat,[0,0,0,(icell-1)*26+2+13,(icell-1)*26+3+13,0.0,socpb])

        for i=0:2
          push!(hopsdat,[0,0,0,(icell-1)*26+5+3*i+2,(icell-1)*26+5+3*i+0+13,-soci,0.0])
          push!(hopsdat,[0,0,0,(icell-1)*26+5+3*i+2,(icell-1)*26+5+3*i+1+13,0.0,soci])
          push!(hopsdat,[0,0,0,(icell-1)*26+5+3*i+0,(icell-1)*26+5+3*i+2+13,soci,0.0])
          push!(hopsdat,[0,0,0,(icell-1)*26+5+3*i+1,(icell-1)*26+5+3*i+2+13,0.0,-soci])
          push!(hopsdat,[0,0,0,(icell-1)*26+5+3*i+0,(icell-1)*26+5+3*i+1,0.0,-soci])
          push!(hopsdat,[0,0,0,(icell-1)*26+5+3*i+0+13,(icell-1)*26+5+3*i+1+13,0.0,soci])
        end

	# extra next nearest neighbor hopping
	if use_pbpb_hop
	push!(hopsdat, [dx,0,0,(icell-1)*26+2,(icellxm-1)*26+2, nnpb,0.0]) #pb x nn hopping
	push!(hopsdat, [0,dy,0,(icell-1)*26+3,(icellym-1)*26+3, nnpb,0.0]) #y
	push!(hopsdat, [0,0,dz,(icell-1)*26+4,(icellzm-1)*26+4, nnpb,0.0]) #z
	push!(hopsdat, [dx,0,0,(icell-1)*26+13+2,(icellxm-1)*26+13+2, nnpb,0.0]) # down spins
	push!(hopsdat, [0,dy,0,(icell-1)*26+13+3,(icellym-1)*26+13+3, nnpb,0.0])  
	push!(hopsdat, [0,0,dz,(icell-1)*26+13+4,(icellzm-1)*26+13+4, nnpb,0.0]) 
	end
	if use_ii_hop1
	push!(hopsdat, [dx,0,0,(icell-1)*26+5, (icellxm-1)*26+5,nni,0.0]) # I nn hopping across Pb
	push!(hopsdat, [0,dy,0,(icell-1)*26+9, (icellym-1)*26+9,nni,0.0])
	push!(hopsdat, [0,0,dz,(icell-1)*26+13, (icellzm-1)*26+13,nni,0.0])
	push!(hopsdat, [dx,0,0,(icell-1)*26+13+5, (icellxm-1)*26+13+5,nni,0.0])
	push!(hopsdat, [0,dy,0,(icell-1)*26+13+9, (icellym-1)*26+13+9,nni,0.0])
	push!(hopsdat, [0,0,dz,(icell-1)*26+13+13, (icellzm-1)*26+13+13,nni,0.0])
	end
	if use_ii_hop2
	push!(hopsdat, [0,dy,0,(icell-1)*26+6, (icellym-1)*26+6,nni2,0.0])#Ix hopping, py to py
	push!(hopsdat, [0,0,dz,(icell-1)*26+7, (icellzm-1)*26+7,nni2,0.0])#Ix hopping, pz to pz
	push!(hopsdat, [dx,0,0,(icell-1)*26+8, (icellxm-1)*26+8,nni2,0.0])#Iy, px to px
	push!(hopsdat, [0,0,dz,(icell-1)*26+10,(icellzm-1)*26+10,nni2,0.0])#Iy, pz to pz
	push!(hopsdat, [dx,0,0,(icell-1)*26+11,(icellxm-1)*26+11,nni2,0.0]) #Iz, px to px 
	push!(hopsdat, [0,dy,0,(icell-1)*26+12,(icellym-1)*26+12,nni2,0.0]) #Iz, py to py
	push!(hopsdat, [0,dy,0,(icell-1)*26+6+13, (icellym-1)*26+6+13,nni2,0.0])#Ix hopping, py to py
	push!(hopsdat, [0,0,dz,(icell-1)*26+7+13, (icellzm-1)*26+7+13,nni2,0.0])#Ix hopping, pz to pz
	push!(hopsdat, [dx,0,0,(icell-1)*26+8+13, (icellxm-1)*26+8+13,nni2,0.0])#Iy, px to px
	push!(hopsdat, [0,0,dz,(icell-1)*26+10+13,(icellzm-1)*26+10+13,nni2,0.0])#Iy, pz to pz
	push!(hopsdat, [dx,0,0,(icell-1)*26+11+13,(icellxm-1)*26+11+13,nni2,0.0]) #Iz, px to px 
	push!(hopsdat, [0,dy,0,(icell-1)*26+12+13,(icellym-1)*26+12+13,nni2,0.0]) #Iz, py to py
	end

      end
    end
  end

  hopsdat = aoa2d(hopsdat)
  #println(hopsdat)
  writedlm("hops.dat",hopsdat)
 
  kpath0 = [ 0.5 0.0 0.0 ; 0.0 0.0 0.0 ; 0.0 0.0 0.5 ]
  nk = 5
  kpath = kpts_path(kpath0,nk)
  kpath = [ 0.5 0.0 0.0 ; 0.4 0.0 0.0 ; 0.3 0.0 0.0 ; 0.2 0.0 0.0; 0.1 0.0 0.0; 0.0 0.0 0.0; 0.0 0.1 0.0 ; 0.0 0.2 0.0; 0.0 0.3 0.0 ; 0.0 0.4 0.0 ; 0.0 0.5 0.0; 0.0 0.0 0.1; 0.0 0.0 0.2; 0.0 0.0 0.3; 0.0 0.0 0.4; 0.0 0.0 0.5]
  nh = ncells*26
  tstr = string(it)
  solvehamk(k->tbhammaker3d(k,nh,hopsdat,spars=true), kpath, [20*ncells-1:20*ncells+2],savewfns=false,verbose=true, spars=true,nev=24,etarget=4.7, outfile="bands$tstr.dat")

  gap = 10.0
  #for ik=1:2*nk
   for ik=1:1
    endat = readdlm(open("work/bands$tstr.dat"), Float64) #readvec(float,"work/k$(ik)/ens")
    #gapt = 0.0
    for ie=1:(size(endat,1)-1)
      if (endat[ie,4]<4.7 && endat[ie+1,4]>4.7)
        gapt= endat[ie+1,4] - endat[ie,4]
        #break
	if (gap> gapt)
	  gap = gapt
	end
      end
    end
    #if (gap> gapt)
     # gap = gapt
    #end
  end
  println(gap) 
  println(Dates.now()) 
  gap #function returns gap
	
end

end 
