subroutine saveout_parameters
  !prints everything and q0 to files
  use global_variables
  use geometries
  use mpi
  use mpi_surf
  implicit none
  
  integer::i,j
  character(len=255) ::filename,ofmt,fname,command
  logical::dir_e
  
  if (myrank==master) then
    fname=trim(outdir)//'/grids/'
    inquire(file=fname,exist=dir_e)    
!     print*, fname
!     print*,dir_e
    if( .not. dir_e) then
      command='mkdir '//trim(fname)
      print*,command
      call system(trim(command))
    end if 

    fname=trim(outdir)//'/grids/gridproperties.dat'
    open(unit=50,file=fname)
    write(unit=50,fmt='("dit ",I5.5)') dit
    write(unit=50,fmt='("nproc ",I3.3)') nproc  
    write(unit=50,fmt='("dims ",I3.3,X,I3.3)') dims(1),dims(2)
    write(unit=50,fmt='("nxi ",I4.4)') nxi
    write(unit=50,fmt='("neta ",I4.4)')neta

    if( (batiopt==0).or.(batiopt==1) )then
      command='cp '//trim(batiname(1))//' '//trim(outdir)//'/.'
      call system(command)
      command='cp '//trim(batiname(2))//' '//trim(outdir)//'/.'
      call system(command)
      command='cp '//trim(batiname(3))//' '//trim(outdir)//'/.'
      call system(command)
      write(unit=50,fmt=*) 'batiname ', trim(batiname(1)),' ', &
			  trim(batiname(2)),' ',trim(batiname(3))
    else if( (batiopt==2).or.(batiopt==3) )then
      command='cp '//trim(batiname(1))//' '//trim(outdir)//'/.'
      call system(command)
      write(unit=50,fmt=*) 'batiname', trim(batiname(1))
    end if
    close(unit=50)       
  end if
  
  write(filename,'("/grids/grid",I3.3,"_",I3.3,".dat")') coords(1),coords(2)
  filename=trim(outdir)//trim(adjustl(filename))
  open(unit=myrank+100,file=filename)
  write(unit=myrank+100,fmt='( I3.3, "  Nbx" )') Nbx
  write(unit=myrank+100,fmt='( I3.3, "  Nby" )') Nby  
  write(unit=myrank+100,fmt='( I3.3, "  si" )') si
  write(unit=myrank+100,fmt='( I3.3, "  ei" )') ei
  write(unit=myrank+100,fmt='( I3.3, "  sj" )') sj
  write(unit=myrank+100,fmt='( I3.3, "  ej" )') ej
  write(unit=myrank+100,fmt='( I3.3, "  coord(1)" )') coords(1)
  write(unit=myrank+100,fmt='( I3.3, "  coord(2)" )') coords(2)
  write(unit=myrank+100,fmt='( I3.3, "  dims(1)" )') dims(1)
  write(unit=myrank+100,fmt='( I3.3, "  dims(2)" )') dims(2)  
  write(unit=myrank+100,fmt='( I3.3, "  rank2d" )') myrank2d
  write(unit=myrank+100,fmt='( I4.3, "  left" )') myleft
  write(unit=myrank+100,fmt='( I4.3, "  right" )') myright
  write(unit=myrank+100,fmt='( I4.3, "  back" )') myback
  write(unit=myrank+100,fmt='( I4.3, "  front" )') myfront
  close(unit=myrank+100)

end subroutine