program Main_Program
    use mtmod
    use choice
    IMPLICIT NONE

! PROGRAM square well simulates a cell
! proliferation along a one dimens-ional
! axis with a square well potential
!
! MODULES mtmod (Mersenne Twister) and
! choice are needed to compile.
!
! Generated : November 2018
! Last updated : 16 June 2020
!
! Created by :
! Celso T. Villano Jr.
! Modified by: 
! Marc Jason P. Bunagan
!========================================

!===========SET VARIABLES=================!
REAL(8)  :: pi, V_pot, &
            bnd_x, r_min, r_max, V_o, step, &
            delta_r, d_crit_grow, d_crit_deform, total_cell_count, &
            choose_cell, decision, r_grow, rand_num, r_ubgrow, &
            rand_num2, r_dist, r_angle, r_d, d_max, d_min
REAL(8) :: start, finish

REAL(8), DIMENSION(2):: r_move
REAL(8), DIMENSION(3):: r_mat
REAL(8), DIMENSION(9) :: kioku   !it means memory
REAL(8), DIMENSION(3000,9) :: A               !memory: 1 row for 1 cell
REAL(8), DIMENSION(100000,2) :: N_vs_t        !t = 100,000

INTEGER :: seed, overlap, overlap_o, &
           A_x, A_y, total_time_step, choose_cell_int, &
           row, counter_cell, counter, overlap_m, m_check, repeated
pi = acos(-1.0)

!===========SET INITIAL PARAMETERS========!
counter = 0 
repeated = 0        !checking progress when running  
overlap = 0
overlap_o = 0
overlap_m = 0
counter_cell = 1    !# of n(t=0) cells

!=====SYSTEM/ALGORITHM PARAMETERS=====!
total_time_step = 100000
total_cell_count = 3000
bnd_x = 10E8       !boundary of environemnt
r_min = 1.0
r_max = (2.0)**(1.0/3.0)*r_min
V_o = vol(r_min)
step = r_min / 6.0
d_max = 2*r_min
delta_r = 0.1*r_min 
d_crit_grow = .5   !r_g = [0,1]
d_crit_deform = .5 !r_d = [0,1]
d_min = 0.03       !Prevent error on r when d is very small

!=======POTENTIALS=========!
!=======SQUARE WELL========!
V_pot = -5          !(epsilon:)

!=====MATRIX (A)=======!
!=======================================================================!
! MATRIX PARAMETERS DEFINITION: 
! Row defines which cell, while column holds the cell's characteristics
! 1st column = x-coordinate 
! 2nd column = y-coordinate
! 3rd column = z-coordinate 
! 4th column = radius
! 5th column = phase of the cell
! 
! If M-phase: these columns are used temporarily for the soon-to-be cells
! 6th column = x-coordinate 
! 7th column = y-coordinate 
! 8th column = z-coordinate
! 9th column = separation distance
!
! 
! 
!=======================================================================!

A_x = total_cell_count
A_y = 9

!+++First Cell+++!
!ALLOCATE(A(A_x,A_y))

A(1,4) = r_min     
A(1,5) = 1.0       !interphase
A(1,1) = 0
A(1,2) = 0
A(1,3) = 0

!+++Allocate Temporary Memory+++!
kioku = 0

!====RANDOMNESS=====!
seed = 1502*xlb_time()
call sgrnd(seed)
!============================================!

!===================MAIN PROGRAM======================!
call cpu_time(start)
do while (counter <= total_time_step)

    counter = counter + 1
    !choose a cell id
    choose_cell = grnd()
    choose_cell_int = xlb_choose(choose_cell, counter_cell)

    !redo the loop if the chosen cell is not initialized
    if (A(choose_cell_int, 4) == 0.0  .AND. A(choose_cell_int, 5) == 0.0) then
        cycle
    end if

!=============================I-PHASE(GROWTH)===========================================!
    if (A(choose_cell_int,5) == 1.0) then                 
        decision = grnd()                                    
        
        if (decision <= d_crit_grow) then                  
            kioku = A(choose_cell_int,:)                    
            rand_num  = grnd()                            
            r_grow    = rand_num*r_min/12.0  !M_cp for cell growth            
            r_ubgrow  = A(choose_cell_int,4) + r_grow
            A(choose_cell_int,4) = increase(A(choose_cell_int,4) ,r_max ,r_grow) 


            do row = 1, counter_cell                    
                if (row == choose_cell_int) then          
                    cycle
                end if
                overlap = overlap_2D(A(choose_cell_int,1), A(choose_cell_int,2), A(choose_cell_int,6), A(choose_cell_int,7), &
                                     A(row,1), A(row,2), A(row,6), A(row,7), &
                                     A(choose_cell_int,4), A(row,4), delta_r)                                                 
                if (overlap == 2) then                                                 
                    A(choose_cell_int,:) = kioku                                      
                    counter = counter - 1
                    exit
                end if
            end do    

            if (overlap < 2) then
                if (A(choose_cell_int,4) >= r_max) then                   
                    A(choose_cell_int,5) =  2.0                        
                    A(choose_cell_int,6) =  A(choose_cell_int,1)          
                    A(choose_cell_int,7) =  A(choose_cell_int,2)  
                    A(choose_cell_int,8) =  A(choose_cell_int,3)    
                    A(choose_cell_int,9) =  d_min                        
                end if
                             
            end if

            
!=============================I-PHASE(MIGRATE)===========================================!
        else                                              
            kioku = A(choose_cell_int,:)               
            rand_num = grnd()                           
            rand_num2 = grnd()                           
            r_dist = rand_num*r_min  !M_cp for cell migration                   
            r_angle = rand_num2 * 2 * pi                
            r_move = migrate(r_dist, r_angle, bnd_x, bnd_x,A(choose_cell_int,1),A(choose_cell_int,2))
            A(choose_cell_int,1) = r_move(1)    
            A(choose_cell_int,2) = r_move(2)                     
         

            do row = 1, counter_cell                  
                if (row == choose_cell_int) then        
                    cycle
                end if
                overlap = overlap_2D(A(choose_cell_int,1), A(choose_cell_int,2), A(choose_cell_int,6), A(choose_cell_int,7), &
                                    A(row,1), A(row,2), A(row,6), A(row,7), &
                                    A(choose_cell_int,4), A(row,4), delta_r)
                if (overlap == 2) then                                           
                    A(choose_cell_int,:) = kioku                                     
                    overlap_m = 2
                    counter = counter - 1   
                    exit
                else if (overlap == 1 .and. overlap_m /= 1) then                                          
                    overlap_m = 1
                end if
            end do 
            
            if (overlap_m < 2) then                            
                do row = 1, counter_cell                    
                    if (row == choose_cell_int) then          
                        cycle
                    end if
                    overlap_o = overlap_2D(kioku(1), kioku(2), kioku(6), kioku(7), &
                                           A(row,1), A(row,2), A(row,6), A(row,7), &
                                           A(choose_cell_int,4), A(row,4), delta_r)
                    if (overlap_o == 1) then                
                        exit
                    end if
                end do 
            end if

            !metroplis algorithm
            rand_num = grnd()
            m_check = metropolis_square(overlap_m,overlap_o,V_pot,rand_num)
            if (m_check == 0) then                                                 
                A(choose_cell_int,:) = kioku                                    
                counter = counter - 1   
            end if
        end if 


!=============================M-PHASE(DEFORM)===========================================!
    else if (A(choose_cell_int,5) == 2.0) then            
        decision = grnd()                                   
        if (decision <= d_crit_deform) then                
            kioku = A(choose_cell_int,:)                   
            rand_num  = grnd()                            
            r_d       = rand_num*.5*r_min  !M_cp for cell deformation                         
            !solve new cell location since sep_d changed                             
            r_mat     = adjust_deform(A(choose_cell_int,1),A(choose_cell_int,6),A(choose_cell_int,9),r_d,d_max)     
            A(choose_cell_int,1) = r_mat(1)                                                          
            A(choose_cell_int,6) = r_mat(2)
            A(choose_cell_int,9) = r_mat(3)                                                        
            A(choose_cell_int,4) = r_deform(A(choose_cell_int,9),V_o)                              
 
            do row = 1, counter_cell                      
                if (row == choose_cell_int) then          
                    cycle
                end if
                overlap = overlap_2D(A(choose_cell_int,1), A(choose_cell_int,2), A(choose_cell_int,6), A(choose_cell_int,7), &
                                    A(row,1), A(row,2), A(row,6), A(row,7), &
                                    A(choose_cell_int,4), A(row,4), delta_r)
                if (overlap == 2) then                                              
                    A(choose_cell_int,:) = kioku                                       
                    counter = counter - 1     
                    exit
                end if
            end do    
            if (overlap < 2) then
                if (A(choose_cell_int,9) >= d_max) then               
                    A(choose_cell_int,9) =  0.0   
                    A(choose_cell_int,5) =  1.0                        
                    A(choose_cell_int,4) =  r_min
                    counter_cell = counter_cell + 1
                    A(counter_cell,1) =  A(choose_cell_int,6)          
                    A(counter_cell,2) =  A(choose_cell_int,7)  
                    A(counter_cell,3) =  A(choose_cell_int,8) 
                    A(counter_cell,4) =  r_min
                    A(counter_cell,5) =  1.0                         
                    A(choose_cell_int,6:9) = 0.0                       
                                                
                end if             
            end if
!=============================M-PHASE(MIGRTATE)===========================================!
        else 
            kioku = A(choose_cell_int,:)             
            rand_num = grnd()                       
            rand_num2 = grnd()                          
            r_dist = rand_num*r_min  !m_cp for cell migration                 
            r_angle = rand_num2 * 2 * pi                
            r_move = migrate(r_dist, r_angle, bnd_x, bnd_x,A(choose_cell_int,1),A(choose_cell_int,2))
            A(choose_cell_int,1) = r_move(1)    
            A(choose_cell_int,2) = r_move(2)
            r_move = migrate(r_dist, r_angle, bnd_x, bnd_x,A(choose_cell_int,6),A(choose_cell_int,7))
            A(choose_cell_int,6) = r_move(1)    
            A(choose_cell_int,7) = r_move(2)                     
         

            do row = 1, counter_cell                    
                if (row == choose_cell_int) then        
                    cycle
                end if
                overlap = overlap_2D(A(choose_cell_int,1), A(choose_cell_int,2), A(choose_cell_int,6), A(choose_cell_int,7), &
                                    A(row,1), A(row,2), A(row,6), A(row,7), &
                                    A(choose_cell_int,4), A(row,4), delta_r)
                if (overlap == 2) then                                                
                    A(choose_cell_int,:) = kioku                                    
                    overlap_m = 2
                    counter = counter - 1   
                    exit
                else if (overlap == 1 .and. overlap_m /= 1) then                                          
                    overlap_m = 1   
                end if
            end do 
            
            if (overlap_m < 2) then                           
                do row = 1, counter_cell                  
                    if (row == choose_cell_int) then         
                        cycle
                    end if
                    overlap_o = overlap_2D(kioku(1), kioku(2), kioku(6), kioku(7), &
                                           A(row,1), A(row,2), A(row,6), A(row,7), &
                                           A(choose_cell_int,4), A(row,4), delta_r)
                    if (overlap_o == 1) then                   
                        exit
                    end if
                end do 
            end if

            !metroplis algorithm
            rand_num = grnd()
            m_check = metropolis_square(overlap_m,overlap_o,V_pot,rand_num)
            if (m_check == 0) then                                               
                A(choose_cell_int,:) = kioku                                    
                counter = counter - 1   
            end if
        end if
    

    else !neither interphase or mitosis? (checking)
        print *, 'cell has wrong data'
    end if


    !checking progress when running
    N_vs_t(counter,:) = (/counter,counter_cell/) 
    if (repeated == 0) then
        repeated = repeated + 5000
    else if (counter == repeated) then
        print *, counter
        repeated = repeated + 5000
    end if

    !Clear some data
    overlap_m = 0
    
end do



!Generate N_vs_t data for Plottong in Python
open(11,file="FileLocation1.txt")
write (11,*) N_vs_t(:,1)

open(12,file="FileLocation2.txt")
write (12,*) N_vs_t(:,2)

!Collect Data of Each Cell
open(21,file="FileLocation3.txt")
write (21,*) A

call cpu_time(finish)
print '("Time = ",f6.3," seconds.")',finish-start
end program Main_Program

!attract 0     = 24-28s
!attract -.29  = 25s
!attract -0.69 = 27s
!attract -1.39 = 32-36s
!attract -5    = 40-44s
!attract -20   = 48s



!deform .1 = 12-14s
!deform .2 = 23-25s
!deform .3 = 28-32s
!deform .4 = 32-38s
!deform  .5 = 40-44s
!deform  .6 = 45-50s
!deform  .7 = 52-56s
!deform  .8 = 67-60s
!deform  .9 = 90-98s

!grow .1 = 12s
!grow .2 = 23-25s
!grow .3 = 32-36s
!grow .4 = 41-42s
!grow .5 = 40-44s
!grow .6 = 45-48s
!grow .7 = 49-51s
!grow .8 = 52-60s
!grow .9 = 76s