module Htools

export solvehamk

#find the bandstructure of a hamiltonian
#hammaker(k) = hamiltonian
#read the wavefunctions using read(wfnfile,Complex128,(dim(hamiltonian),nbands))
#gives evecs[ig,ibands]
# Options:
# savewfns: whether to save the wavefunctions
# wfnfmt : print in binary (false) or ascii (true)
# dir: output directory
# spars: use sparse matrix diagonalization, 
# nev: number of eigenvalues to find in sparse matrix diagonalization
# etarget: find eigenvalues close to this in sparse matrix diagonalization
#
#
# if(wfnfmt), the wfn file has the format
# re(band1,g1) im(band1,g1) re(band2,g1) im(band2,g1) ...
# re(band1,g2) im(band1,g2) re(band2,g2) im(band2,g2) ...
# ...
#
# Additionally, prints to file
# k.dat: kpoints in invang
# bands.dat: kpoints and energies
function solvehamk(hammaker, kpts, bands; savewfns=true,wfnfmt=false,verbose=true, dir="work", spars=false,nev=8,etarget=1.0,outfile="bands.dat")

  run(`mkdir -p $dir`)
  ndim = size(kpts,2)
  
  fk = open("$dir/k.dat","w")
  for ik = 1:size(kpts,1)
    kpt = vec(kpts[ik,:])
    for idim = 1:ndim
      @printf(fk,"%s ",kpt[idim])
    end
    @printf(fk,"\n")
  end
  close(fk)

  fout = open("$dir/$outfile","w")
  for ik = 1:size(kpts,1)
    if(verbose)
      println("doing ik: ",ik)
    end
    kpt = vec(kpts[ik,:])
    if(savewfns)
      run(`mkdir -p $dir/k$ik`)
    end
    if(spars)
      ens,evecs = eigs(hammaker(kpt); nev=nev, sigma=etarget)
      #ens = real(ens)
      bands = [1:nev]
    else
	#print(hammaker(kpt))
	ens, evecs = eigs(hammaker(kpt))
    end
    perm = sortperm(ens,by=real)
    if(savewfns)
      writedlm("$dir/k$ik/ens", ens[perm][bands], ' ' ) 
      evecs_sorted = evecs[:,perm]
      if(wfnfmt)
        nh = size(evecs_sorted,1)
        nb = size(bands,1)
        evecs_towrite = zeros(typeof(real(evecs_sorted[1,1])),(nh,2*nb))
        for ig=1:nh
          for ib=1:nb
            ib1 = bands[ib]
            evecs_towrite[ig,2*(ib-1)+1]=real(evecs_sorted[ig,ib1])
            evecs_towrite[ig,2*(ib-1)+2]=imag(evecs_sorted[ig,ib1])
          end
        end
        writedlm("$dir/k$ik/evecs", evecs_towrite, ' ' )
      else
        wfnfile = open("$dir/k$ik/evecs","w")
        write(wfnfile, evecs_sorted[:,bands])
        close(wfnfile)
      end
    end
    for b in sort(real(ens))
      for idim = 1:ndim
        @printf(fout,"%s ",kpt[idim])
      end
      @printf(fout,"%s\n",b)
    end
  end
  close(fout)
end

end
