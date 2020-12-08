module choice
    use mtmod
    use iso_fortran_env, only: r128 => real128
    !==============================================
    ! MODULE choice holds the list of functions
    ! that supplements the programs harmonic well
    ! and square well.
    !
    ! Generated : November 2018
    ! Last edited : April 2020
    !
    ! Created by :
    ! Celso T. Villano Jr.
    ! Modified for 2D by:
    ! Marc Jason P. Bunagan
    !
    !==============================================
    contains
    function xlb_time() result(t) !#nochange
        integer, dimension (8) :: values
        call date_and_time (VALUES=values)
        t = values(8)
    end function

    !FUNCTION: Increase radius by rand_inc #nochange
    function increase(x ,r_max ,rand_inc) result (z)
        real(8), intent(in) :: rand_inc, r_max, x
        real(8) :: new, z
        new = x + rand_inc
        if (new < r_max ) then
            z = new
        else
            z = r_max
        end if
    end function

    !FUNCTION: Decrease radius by rand_inc #nochange
    function decrease(x , r_ubgrow) result (z)
        real(8), intent(in) :: x , r_ubgrow
        real(8) :: z
        z = x - r_ubgrow
    end function

    !Migrate cell until it reaches the boundary #changeto2D
    function migrate(r_dist, r_angle, bnd_x, bnd_y,x,y) result(r)
        real(8), intent(in) :: r_dist, r_angle, bnd_x, bnd_y, x, y
        real(8) :: rx, ry, x_coord, y_coord
        real(8), dimension(2):: r
        x_coord = r_dist * COS(r_angle)
        y_coord = r_dist * SIN(r_angle) 
        rx = 0.0
        ry = 0.0
        if (x + x_coord < bnd_x .AND. x + x_coord > -1.0*bnd_x  .AND. y + y_coord < bnd_y .AND. y + y_coord > -1.0*bnd_y) then
            rx = x + x_coord
            ry = y + y_coord
        else
            rx = x_coord
            ry = y_coord
            print *, "cell is at boundary"
        end if
        r(1) = rx 
        r(2) = ry
    end function

    !Deform the dividing cells in such a way to preserve its volume #precisionchanged
    function r_deform(dt2,V_o2) result(rt)
        real(8), intent(in) :: dt2, V_o2

        real(r128) :: dt, V_o
        real(r128) :: pi, qt1, qt2, qt, rt2
        real(r128) :: ut 
        real(8) :: rt

        dt  = dt2
        V_o = V_o2
        pi  = acos(-1.0_r128)
        qt1 = (-1.0_r128)*pi*(dt**(3.0_r128))-(48.0_r128)*V_o 
        qt2 = (4.0_r128*SQRT(6.0_r128)*SQRT(24.0_r128*V_o**2.0_r128+(dt**(3.0_r128))*pi*V_o))
        qt  = qt1 + qt2     
        qt  = abs(qt)
        ut  =  (-1.0_r128*pi**(1.0_r128/3.0_r128))*((1.0_r128/qt)**(1.0_r128/3.0_r128))
        rt2  = (-1.0_r128/4.0_r128)* (dt+dt**2.0_r128*ut+1.0_r128/ut)
        rt =   rt2
    end function
    
    !Vol_sphere #nochange
    function vol(r) result(V)
        real(8) , intent(in) :: r 
        real(8) :: V, pi 
        pi = ACOS(-1.0)
        V = (4.0/3.0)*pi*r**3.0
    end function

    !Area_sphere #new
    function area(r) result(V)
        real(8) , intent(in) :: r 
        real(8) :: V, pi 
        pi = ACOS(-1.0)
        V = pi*r**2.0
    end function

    !new code
    function overlap_2D(r_cix, r_ciy, r_cmx, r_cmy, &
                        r_rix, r_riy, r_rmx, r_rmy, &
                        r_cl, r_rl, delta_r) result(overlap)
        !r_c -> chosen cell, r_r -> do row cell
        real(8), intent(in) :: r_cix, r_ciy, r_cmx, r_cmy
        real(8), intent(in) :: r_rix, r_riy, r_rmx, r_rmy
        real(8), intent(in) :: r_cl, r_rl, delta_r  
        real(8) :: r_small
        real(8), dimension(4) :: r_cont
        integer :: overlap
        !overlap -> 0 = no overlap, 2 = hard core overlap, 1 = membrane overlap
        !i-cell -> i-cell [r_cont(1)]
        r_cont(1) = ((r_cix - r_rix)**2 + (r_ciy - r_riy)**2)**(1./2.) 
        !i-cell -> m-cell [r_cont(2)]
        if  (r_cmx /= 0) then
            r_cont(2) = ((r_cmx - r_rix)**2 + (r_cmy - r_riy)**2)**(1./2.)
        else
            r_cont(2) = r_cont(1) + 1
        end if
        !m-cell -> i-cell [r_cont(3)]
        if  (r_rmx /= 0) then
            r_cont(3) = ((r_cix - r_rmx)**2 + (r_ciy - r_rmy)**2)**(1./2.)
        else
            r_cont(3) = r_cont(1) + 1
        end if
        !m-cell -> m-cell [r_cont(4)]
        if  (r_cmx /= 0 .and. r_rmx /= 0) then
            r_cont(4) = ((r_cmx - r_rmx)**2 + (r_cmy - r_rmy)**2)**(1./2.)
        else
            r_cont(4) = r_cont(1) + 1
        end if
        !choose the smallest distance
        r_small = minval(r_cont)
        !decide the kind of overlap
        if (r_small > r_cl + r_rl + 2*delta_r) then
            overlap = 0
        else if (r_small > r_cl + r_rl) then
            overlap = 1
        else 
            overlap = 2
        end if
    end function


    !used to choose a cell #nochange
    function xlb_choose(rand_num, total_cell_count_int) result(choose_cell_int)
        real(8), intent(in) :: rand_num
        integer, intent(in) :: total_cell_count_int
        real(8) :: a, total_cell_count
        INTEGER :: choose_cell_int
        total_cell_count = total_cell_count_int
        a = rand_num * total_cell_count
        choose_cell_int = INT(a) + 1
    end function


    !used to adjust mitotic cell deformation in 1D #new
    function adjust_deform(x1,x2,d_old, d_new, d_max) result (r_mat)
        real(8), intent(in) :: x1, x2, d_new, d_old, d_max
        real(8) :: d_now
        real(8), dimension(3) :: r_mat
        d_now = d_new + d_old
        if (d_now < d_max) then
            r_mat(1) = x1 - d_new/2
            r_mat(2) = x2 + d_new/2
            r_mat(3) = d_now
        else
            r_mat(1) = x1 + d_old/2 - d_max/2
            r_mat(2) = x2 - d_old/2 + d_max/2
            r_mat(3) = d_max
        end if
    end function

    !Metropolis for a Square Potential #new  !1 = to accept, 0 = to reject
    function metropolis_square(overlap_m, overlap_o, V_pot, rand_num) result(m_check)
        integer, intent(in) :: overlap_m, overlap_o
        real(8), intent(in) :: V_pot, rand_num
        integer :: m_check
        real(8) :: e
        e = EXP(1.0)

        if (overlap_m == 0 .and. overlap_o == 1) then
            if (rand_num < e**(V_pot)) then
                m_check = 1
            else
                m_check = 0
            end if
        else
            m_check = 1
        end if

    end function
end module

    
        






