module max_parameters
  
  implicit none
  
  integer, parameter :: maxatoms = 15000, maxCoo = 3*maxatoms, &
       maxfpi = 1200, maxfpicoo = 3*maxfpi, maxcomp = 10, &
       maxnuptdl = 10000, maxtgr = 21, maxpotch = 50, &
       maxprs = 150000, maxprscoo = 3*maxprs
  
end module max_parameters
