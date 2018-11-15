subroutine residual_func(snes_in, stateVector, residualVector, ctx, localerr)
#include <petsc/finclude/petscsnes.h>
   use petscsnes
   use contexts
   use grids
   use input
   use geometry
   use mp
   use diffusion
   use source
   implicit none
   SNES:: snes_in
   PetscScalar,dimension(:),intent(in):: stateVector
   PetscScalar,dimension(:),intent(inout):: residualVector
   type(resContext):: ctx
   PetscErrorCode:: localerr
   real,dimension(:),allocatable:: Efield
   integer :: idx, ir, ip, ix, idx_rp1, idx_rm1, idx_pp1, idx_pm1, idx_xp1, idx_xm1
   real :: delta_t, t, temp
   real :: jacob_rl,jacob_rr,jacob_pl,jacob_pr,jacob_xl,jacob_xr
   real :: jacob_rm1,jacob_rp1,jacob_pm1,jacob_pp1,jacob_xm1,jacob_xp1, jacob_c
   real :: Dr_l, Dr_r, Dp_l, Dp_r, Dx_l, Dx_r, Ar_l, Ar_r, Ap_l, Ap_r, Ax_l, Ax_r

   allocate(Efield(Nr))

   select case (efield_option)
   case ("none")
      Efield(:) = 0.0
   case ("inductive")
      ! Get electric field from its home processor
      ! TODO
   end select

   t = ctx%time
   delta_t = ctx%delta_t

   do idx = firstLocalRow, lastLocalRow
      ir = get_idx_r(idx)
      ip = get_idx_p(idx)
      ix = get_idx_x(idx)

      idx_rp1 = get_idx(ir+1,ip,ix)
      idx_rm1 = get_idx(ir-1,ip,ix)
      idx_pp1 = get_idx(ir,ip+1,ix)
      idx_pm1 = get_idx(ir,ip-1,ix)
      idx_xp1 = get_idx(ir,ip,ix+1)
      idx_xm1 = get_idx(ir,ip,ix-1)

      ! Coordinate Jacobian defined in cell center
      jacob_c = jacob(rgrid(ir),pgrid(ip),xgrid(ix))

      temp = 0.0

      if (ir == 1) then
         jacob_rr = jacob(rgrid_edge(ir+1),pgrid(ip),xgrid(ix))
         jacob_rp1 = jacob(rgrid(ir+1),pgrid(ip),xgrid(ix))

         Ar_r = jacob_rr*Ar(rgrid_edge(ir+1),pgrid(ip),xgrid(ix),t)
         Dr_r = jacob_c*Drr(rgrid(ir),pgrid(ip),pgrid(ix),t)*jacob_rp1*Drr(rgrid(ir+1),pgrid(ip),xgrid(ix),t) * ( rgrid(ir+1) - rgrid(ir) ) / & 
            ( (rgrid(ir+1)-rgrid_edge(ir+1))*jacob_c*Drr(rgrid(ir),pgrid(ip),xgrid(ix),t) + (rgrid_edge(ir+1)-rgrid(ir))*jacob_rp1*Drr(rgrid(ir+1),pgrid(ip),xgrid(ix),t) )

         temp=temp+ Ar_r * stateVector( get_idx(ir_upwind(ir+1,ip,ix),ip,ix) ) / (rgrid_edge(ir+1)-rgrid_edge(ir))
         temp=temp- (Dr_r/( (rgrid(ir+1)-rgrid(ir))*(rgrid_edge(ir+1)-rgrid_edge(ir)))) * (stateVector(idx_rp1) - stateVector(idx))
      else if (ir == Nr) then
         jacob_rl = jacob(rgrid_edge(ir),pgrid(ip),xgrid(ix))
         jacob_rr = jacob(rgrid_edge(ir+1),pgrid(ip),xgrid(ix))
         jacob_rm1 = jacob(rgrid(ir-1),pgrid(ip),xgrid(ix))

         Ar_l = jacob_rl*Ar(rgrid_edge(ir),pgrid(ip),xgrid(ix),t)
         Ar_r = jacob_rr*Ar(rgrid_edge(ir+1),pgrid(ip),xgrid(ix),t)

         Dr_l = jacob_rm1*Drr(rgrid(ir-1),pgrid(ip),pgrid(ix),t)*jacob_c*Drr(rgrid(ir),pgrid(ip),xgrid(ix),t) * ( rgrid(ir) - rgrid(ir-1) ) / & 
            ( (rgrid(ir)-rgrid_edge(ir)) * jacob_rm1*Drr(rgrid(ir-1),pgrid(ip),xgrid(ix),t) + (rgrid_edge(ir+1)-rgrid(ir))*jacob_c*Drr(rgrid(ir),pgrid(ip),xgrid(ix),t) )

         Dr_r = jacob_c*Drr(rgrid(ir),pgrid(ip),pgrid(ix),t)*jacob_rp1*Drr(rgrid(ir+1),pgrid(ip),xgrid(ix),t) * ( rgrid(ir+1) - rgrid(ir) ) / & 
            ( (rgrid(ir+1)-rgrid_edge(ir+1))*jacob_c*Drr(rgrid(ir),pgrid(ip),xgrid(ix),t) + (rgrid_edge(ir+1)-rgrid(ir))*jacob_rp1*Drr(rgrid(ir+1),pgrid(ip),xgrid(ix),t) )

         temp=temp- Ar_l * stateVector( get_idx(ir_upwind(ir,ip,ix),ip,ix) ) / (rgrid_edge(ir+1)-rgrid_edge(ir))
         temp=temp+ 0.0

         ! Use a "ghost" cell with vanishing state vector. This cell has a width equal to the last r cell
         temp=temp+ (Dr_l/( (rgrid(ir)-rgrid(ir-1))*(rgrid_edge(ir+1)-rgrid_edge(ir)))) * ( stateVector(idx) - stateVector(idx_rm1) )
         temp=temp- (Dr_r/( (rgrid(ir)-rgrid(ir-1))*(rgrid_edge(ir+1)-rgrid_edge(ir)))) * (0.0 - stateVector(idx))

      else
         jacob_rl = jacob(rgrid_edge(ir),pgrid(ip),xgrid(ix))
         jacob_rr = jacob(rgrid_edge(ir+1),pgrid(ip),xgrid(ix))
         jacob_rm1 = jacob(rgrid(ir-1),pgrid(ip),xgrid(ix))
         jacob_rp1 = jacob(rgrid(ir+1),pgrid(ip),xgrid(ix))

         Ar_l = jacob_rl*Ar(rgrid_edge(ir),pgrid(ip),xgrid(ix),t)
         Ar_r = jacob_rr*Ar(rgrid_edge(ir+1),pgrid(ip),xgrid(ix),t)

         Dr_l = jacob_rm1*Drr(rgrid(ir-1),pgrid(ip),pgrid(ix),t)*jacob_c*Drr(rgrid(ir),pgrid(ip),xgrid(ix),t) * ( rgrid(ir) - rgrid(ir-1) ) / & 
            ( (rgrid(ir)-rgrid_edge(ir)) * jacob_rm1*Drr(rgrid(ir-1),pgrid(ip),xgrid(ix),t) + (rgrid_edge(ir+1)-rgrid(ir))*jacob_c*Drr(rgrid(ir),pgrid(ip),xgrid(ix),t) )

         Dr_r = jacob_c*Drr(rgrid(ir),pgrid(ip),pgrid(ix),t)*jacob_rp1*Drr(rgrid(ir),pgrid(ip),xgrid(ix),t) * ( rgrid(ir) - rgrid(ir-1) ) / & 
            ( (rgrid(ir)-rgrid_edge(ir))*jacob_c*Drr(rgrid(ir),pgrid(ip),xgrid(ix),t) + (rgrid_edge(ir+1)-rgrid(ir))*jacob_rp1*Drr(rgrid(ir),pgrid(ip),xgrid(ix),t) )

         temp=temp- Ar_l * stateVector( get_idx(ir_upwind(ir,ip,ix),ip,ix) ) / (rgrid_edge(ir+1)-rgrid_edge(ir))
         if (ir_upwind(ir+1,ip,ix) > Nr) then
            temp = temp+0.0
         else 
            temp=temp+ Ar_r * stateVector( get_idx(ir_upwind(ir+1,ip,ix),ip,ix) ) / (rgrid_edge(ir+1)-rgrid_edge(ir))
         end if

         temp=temp+ (Dr_l/( (rgrid(ir)-rgrid(ir-1))*(rgrid_edge(ir+1)-rgrid_edge(ir)))) * ( stateVector(idx) - stateVector(idx_rm1) )
         temp=temp- (Dr_r/( (rgrid(ir+1)-rgrid(ir))*(rgrid_edge(ir+1)-rgrid_edge(ir)))) * (0.0 - stateVector(idx))
      end if

      if (ip == 1) then
         jacob_pr = jacob(rgrid(ir),pgrid_edge(ip+1),xgrid(ix))
         jacob_pp1 = jacob(rgrid(ir),pgrid(ip+1),xgrid(ix))

         Ap_r = jacob_pr*Ar(rgrid(ir),pgrid_edge(ip+1),xgrid(ix),t) 
         Dp_r = jacob_c*Dpp(rgrid(ir),pgrid(ip),pgrid(ix),t)*jacob_pp1*Dpp(rgrid(ir),pgrid(ip+1),xgrid(ix),t) * ( pgrid(ip+1) - pgrid(ip) ) / & 
            ( (pgrid(ip+1)-pgrid_edge(ip+1))*jacob_c*Dpp(rgrid(ir),pgrid(ip),xgrid(ix),t) + (pgrid_edge(ip+1)-pgrid(ip))*jacob_pp1*Dpp(rgrid(ir),pgrid(ip+1),xgrid(ix),t) )

         if (ip_upwind(ir,ip+1,ix) < 1) then
            temp=temp+ 0.0
         else
            temp=temp+ Ap_r * stateVector( get_idx(ir,ip_upwind(ir,ip+1,ix),ix) ) / (pgrid_edge(ip+1)-pgrid_edge(ip))
         end if
         ! Use a "ghost" cell with vanishing state vector. This cell has a width equal to the last p cell
         temp=temp- (Dp_r/( (pgrid(ip)-pgrid(ip-1))*(pgrid_edge(ip+1)-pgrid_edge(ip)))) * (0.0 - stateVector(idx))
      else if (ip == Np) then
         jacob_pl = jacob(rgrid(ir),pgrid_edge(ip),xgrid(ix))
         jacob_pr = jacob(rgrid(ir),pgrid_edge(ip+1),xgrid(ix))
         jacob_pm1 = jacob(rgrid(ir),pgrid(ip-1),xgrid(ix))
         jacob_pp1 = jacob(rgrid(ir),pgrid(ip+1),xgrid(ix))

         Ap_l = jacob_pl*Ar(rgrid(ir),pgrid_edge(ip),xgrid(ix),t) 
         Ap_r = jacob_pr*Ar(rgrid(ir),pgrid_edge(ip+1),xgrid(ix),t) 

         Dp_l = jacob_pm1*Dpp(rgrid(ir),pgrid(ip-1),pgrid(ix),t)*jacob_c*Dpp(rgrid(ir),pgrid(ip),xgrid(ix),t) * ( pgrid(ip) - pgrid(ip-1) ) / & 
            ( (pgrid(ip)-pgrid_edge(ip))*jacob_pm1*Dpp(rgrid(ir),pgrid(ip-1),xgrid(ix),t) + (pgrid_edge(ip+1)-pgrid(ip))*jacob_c**Dpp(rgrid(ir),pgrid(ip),xgrid(ix),t) )

         Dp_r = jacob_c*Dpp(rgrid(ir),pgrid(ip),pgrid(ix),t)*jacob_pp1*Dpp(rgrid(ir),pgrid(ip),xgrid(ix),t) * ( pgrid(ip) - pgrid(ip-1) ) / & 
            ( (pgrid(ip)-pgrid_edge(ip))*jacob_c*Dpp(rgrid(ir),pgrid(ip),xgrid(ix),t) + (pgrid_edge(ip+1)-pgrid(ip))*jacob_pp1*Dpp(rgrid(ir),pgrid(ip),xgrid(ix),t) )

         temp=temp- Ap_l * stateVector( get_idx(ir,ip_upwind(ir,ip,ix),ix) ) / (pgrid_edge(ip+1)-pgrid_edge(ip))
         if (ip_upwind(ir,ip+1,ix) > Np) then
            temp=temp+0.0
         else
            temp=temp+ Ap_r * stateVector( get_idx(ir,ip_upwind(ir,ip+1,ix),ix) ) / (pgrid_edge(ip+1)-pgrid_edge(ip))
         end if

         temp=temp+ (Dp_l/( (pgrid(ip)-pgrid(ip-1))*(pgrid_edge(ip+1)-pgrid_edge(ip)))) * ( stateVector(idx) - stateVector(idx_pm1) )
         temp=temp- (Dp_r/( (pgrid(ip+1)-pgrid(ip))*(pgrid_edge(ip+1)-pgrid_edge(ip)))) * (0.0 - stateVector(idx))

      else
         jacob_pl = jacob(rgrid(ir),pgrid_edge(ip),xgrid(ix))
         jacob_pr = jacob(rgrid(ir),pgrid_edge(ip+1),xgrid(ix))
         jacob_pm1 = jacob(rgrid(ir),pgrid(ip-1),xgrid(ix))
         jacob_pp1 = jacob(rgrid(ir),pgrid(ip+1),xgrid(ix))

         Ap_l = jacob_pl*Ar(rgrid(ir),pgrid_edge(ip),xgrid(ix),t) 
         Ap_r = jacob_pr*Ar(rgrid(ir),pgrid_edge(ip+1),xgrid(ix),t) 

         Dp_l = jacob_pm1*Dpp(rgrid(ir),pgrid(ip-1),pgrid(ix),t)*jacob_c*Dpp(rgrid(ir),pgrid(ip),xgrid(ix),t) * ( pgrid(ip) - pgrid(ip-1) ) / & 
            ( (pgrid(ip)-pgrid_edge(ip))*jacob_pm1*Dpp(rgrid(ir),pgrid(ip-1),xgrid(ix),t) + (pgrid_edge(ip+1)-pgrid(ip))*jacob_c**Dpp(rgrid(ir),pgrid(ip),xgrid(ix),t) )

         Dp_r = jacob_c*Dpp(rgrid(ir),pgrid(ip),pgrid(ix),t)*jacob_pp1*Dpp(rgrid(ir),pgrid(ip+1),xgrid(ix),t) * ( pgrid(ip+1) - pgrid(ip) ) / & 
            ( (pgrid(ip+1)-pgrid_edge(ip+1))*jacob_c*Dpp(rgrid(ir),pgrid(ip),xgrid(ix),t) + (pgrid_edge(ip+1)-pgrid(ip))*jacob_pp1*Dpp(rgrid(ir),pgrid(ip+1),xgrid(ix),t) )
         temp=temp- Ap_l * stateVector( get_idx(ir,ip_upwind(ir,ip,ix),ix) ) / (pgrid_edge(ip+1)-pgrid_edge(ip))
         temp=temp+ Ap_r * stateVector( get_idx(ir,ip_upwind(ir,ip+1,ix),ix) ) / (pgrid_edge(ip+1)-pgrid_edge(ip))

         temp=temp+ (Dp_l/( (pgrid(ip)-pgrid(ip-1))*(pgrid_edge(ip+1)-pgrid_edge(ip)))) * ( stateVector(idx) - stateVector(idx_pm1) )
         temp=temp- (Dp_r/( (pgrid(ip+1)-pgrid(ip))*(pgrid_edge(ip+1)-pgrid_edge(ip)))) * (stateVector(idx_pp1) - stateVector(idx))
      end if

      if (ix == 1) then
         jacob_xr = jacob(rgrid(ir),pgrid(ip),xgrid_edge((ix+1)))
         jacob_xp1 = jacob(rgrid(ir),pgrid(ip),xgrid(ix+1))

         Ax_r = jacob_xr*Ar(rgrid(ir),pgrid(ip),xgrid_edge(ix+1),t) 
         Dx_r = jacob_c*Dxx(rgrid(ir),pgrid(ip),pgrid(ix),t)*jacob_xp1*Dxx(rgrid(ir),pgrid(ip),xgrid(ix+1),t) * ( xgrid(ix+1) - xgrid(ix) ) / & 
            ( (xgrid(ix+1)-xgrid_edge(ix+1))*jacob_c*Dxx(rgrid(ir),pgrid(ip),xgrid(ix),t) + (xgrid_edge(ix+1)-xgrid(ix))*jacob_xp1*Dxx(rgrid(ir),pgrid(ip),xgrid(ix+1),t) )

         temp=temp+ Ax_r * stateVector( get_idx(ir,ip,ix_upwind(ir,ip,ix+1)) ) / (xgrid_edge(ix+1)-xgrid_edge(ix))
         temp=temp- (Dx_r/( (xgrid(ix+1)-xgrid(ix))*(xgrid_edge(ix+1)-pgrid_edge(ix)))) * (stateVector(idx_xp1) - stateVector(idx))
      else if (ix == Nx) then
         jacob_xl = jacob(rgrid(ir),pgrid(ip),xgrid_edge(ix))
         jacob_xm1 = jacob(rgrid(ir),pgrid(ip),xgrid(ix-1))

         Ax_l = jacob_xl*Ar(rgrid(ir),pgrid(ip),xgrid_edge(ix),t) 
         Dx_l = jacob_xm1*Dxx(rgrid(ir),pgrid(ip),pgrid(ix-1),t)*jacob_c*Dxx(rgrid(ir),pgrid(ip),xgrid(ix),t) * ( xgrid(ix) - xgrid(ix-1) ) / & 
            ( (xgrid(ix)-xgrid_edge(ix))*jacob_xm1*Dxx(rgrid(ir),pgrid(ip),xgrid(ix-1),t) + (xgrid_edge(ix+1)-xgrid(ix))*jacob_c*Dxx(rgrid(ir),pgrid(ip),xgrid(ix),t) )

         temp=temp- Ax_l * stateVector( get_idx(ir,ip,ix_upwind(ir,ip,ix)) ) / (xgrid_edge(ix+1)-xgrid_edge(ix))
         temp=temp+ (Dx_l/( (xgrid(ix)-xgrid(ix-1))*(xgrid_edge(ix+1)-pgrid_edge(ix)))) * ( stateVector(idx) - stateVector(idx_xm1) )
      else
         jacob_xl = jacob(rgrid(ir),pgrid(ip),xgrid_edge(ix))
         jacob_xr = jacob(rgrid(ir),pgrid(ip),xgrid_edge((ix+1)))

         jacob_xm1 = jacob(rgrid(ir),pgrid(ip),xgrid(ix-1))
         jacob_xp1 = jacob(rgrid(ir),pgrid(ip),xgrid(ix+1))

         ! Pre-defining the advection and diffusion coefficients. Will probably save this to arrays later
         Ax_l = jacob_xl*Ar(rgrid(ir),pgrid(ip),xgrid_edge(ix),t) 
         Ax_r = jacob_xr*Ar(rgrid(ir),pgrid(ip),xgrid_edge(ix+1),t) 

         Dx_l = jacob_xm1*Dxx(rgrid(ir),pgrid(ip),pgrid(ix-1),t)*jacob_c*Dxx(rgrid(ir),pgrid(ip),xgrid(ix),t) * ( xgrid(ix) - xgrid(ix-1) ) / & 
            ( (xgrid(ix)-xgrid_edge(ix))*jacob_xm1*Dxx(rgrid(ir),pgrid(ip),xgrid(ix-1),t) + (xgrid_edge(ix+1)-xgrid(ix))*jacob_c*Dxx(rgrid(ir),pgrid(ip),xgrid(ix),t) )

         Dx_r = jacob_c*Dxx(rgrid(ir),pgrid(ip),pgrid(ix),t)*jacob_xp1*Dxx(rgrid(ir),pgrid(ip),xgrid(ix+1),t) * ( xgrid(ix+1) - xgrid(ix) ) / & 
            ( (xgrid(ix+1)-xgrid_edge(ix+1))*jacob_c*Dxx(rgrid(ir),pgrid(ip),xgrid(ix),t) + (xgrid_edge(ix+1)-xgrid(ix))*jacob_xp1*Dxx(rgrid(ir),pgrid(ip),xgrid(ix+1),t) )

         temp=temp- Ax_l * stateVector( get_idx(ir,ip,ix_upwind(ir,ip,ix)) ) / (xgrid_edge(ix+1)-xgrid_edge(ix))
         temp=temp+ Ax_r * stateVector( get_idx(ir,ip,ix_upwind(ir,ip,ix+1)) ) / (xgrid_edge(ix+1)-xgrid_edge(ix))

         temp=temp+ (Dx_l/( (xgrid(ix)-xgrid(ix-1))*(xgrid_edge(ix+1)-pgrid_edge(ix)))) * ( stateVector(idx) - stateVector(idx_xm1) )
         temp=temp- (Dx_r/( (xgrid(ix+1)-xgrid(ix))*(xgrid_edge(ix+1)-pgrid_edge(ix)))) * (stateVector(idx_xp1) - stateVector(idx))
      end if

      ! Coordinate jacobian defined on cell faces

      !      residualVector(idx) = temp + ( (stateVector(idx) - ctx%prevStateVector(idx))/delta_t) - sourcefunc(rgrid(ir),pgrid(ip),xgrid(ix),t)
            residualVector(idx) = temp + ( (stateVector(idx) - 0.0) - sourcefunc(rgrid(ir),pgrid(ip),xgrid(ix),t))
   end do
   
end subroutine residual_func

subroutine jacobian_func(snes_in, stateVector, matrix_local, matrix_pc, context, localerr)
#include <petsc/finclude/petscsnes.h>
   use petscsnes
   use contexts
   implicit none
   SNES:: snes_in
   Vec,intent(in):: stateVector
   Mat:: matrix_local, matrix_pc
   type(resContext):: context
   PetscErrorCode:: localerr

end subroutine jacobian_func

