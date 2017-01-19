				!Get v dot E
				vdotE_inst = 0.
				call djpost_jf(part,vdotE_inst,npp,noff,qme,0.,tdjpost,nx,ny,1,inorder,djopt)
				! add guard cells for current in x direction
				call aguard(vdotE_inst,nyp,nx,inorder)
				! add guard cells for current in y direction
				call paguard(vdotE_inst,kstrt,nvp,nx,kyp,ngds)	
				!vdotE_inst now holds vE(t-dt)
				vdotE_sum = vdotE_sum + vdotE_inst
!---
				! Get u. In differential form, du/dt = ( u(t)-u(t-2dt) )/2dt, this is centered at t-dt
				!but when summed for integral, then only need u(t) and u(t-dt) because all other terms in
				!sum cancel.  So, when this diag is dumped, we are dumping: 0.5*( u(t)+u(t-dt) )
				!Total integral from t=0 to t-dt is 0.5*( u(t)+u(d-dt) - (u(0)+u(dt)) )
				!!!NOTE that this looks like it's at t, but it's really the integral from t=0 to t-dt without
				!the initial conditions!!!  Sum the terms of the differenced derivative written above and you'll see.
				!use grad phi because fxye has unwanted smoothing in it.
				grad_phi_dt2 = grad_phi_dt
				grad_phi_dt = grad_phi		!Store old grad_phi
	
				!Find field from gradient of potential
				! calculate potential in fourier space
				call pois(qt,sfieldt,ffc,temp_we,nx,ny,kstrt)
				call ipgradf2(sfieldt,grad_phit,nx,ny,kstrt)
				!Invert fft to get grad_phi in real space: grad_phi(t)
				isign = 1
				call fft(grad_phi,grad_phit,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
				! copy guard cells for current in y direction
				call pcguard(grad_phi,kstrt,nvp,kyp,inorder)
				! copy guard cells for current in x direction
				call cguard(grad_phi,nyp,nx,inorder)
				
				!grad_phi = grad_phi * (-1.)  This step is unneccesary cause all we want is grad_phi**2
				
				!grad_phi_dt = 0.5*( (Ex(t)**2+Ex(t-dt)**2)/2. + (Ey(t)**2+Ey(t-dt)**2)/2. )




!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ NEED THIS FOR INTEGRAL! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!				grad_phi_dt =  0.25*(grad_phi**2 + grad_phi_dt**2)





!---
				!Now get ES Poynting vector
				!Use Viktor's current deposit
				! the 0. is dth=0. so it doesn't do any push
				! the 1 is ipbc=1, which means periodic 2d boundary condition
				cu = 0.


				call djpost(part,cu,npp,noff,qme,0.,tdjpost,nx,ny,1,inorder,djopt)
				! add guard cells for current in x direction
				call aguard(cu,nyp,nx,inorder)
				! add guard cells for current in y direction
				call paguard(cu,kstrt,nvp,nx,kyp,ngds)
				! transform current to fourier space
				isign = -1
				call fft(cu,cut,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
				! smooth current
				call sbpois(cut,tempt,ffc,nx,ny,kstrt)
				! transform current to back
				isign = 1
				call fft(cu,tempt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
				! copy guard cells for current in y direction
				call pcguard(cu,kstrt,nvp,kyp,inorder)
				! copy guard cells for current in x direction
				call cguard(cu,nyp,nx,inorder)
				
				cuold = (cu + cuold)*0.5	!cuold gets j(t-dt) = ( j(t-dt/2) + j(t-3dt/2) )/2

!				cuold = 0.
!				call djpost_jf(part,cuold,npp,noff,qme,0.,tdjpost,nx,ny,1,inorder,djopt)
				! add guard cells for current in x direction
!				call aguard(cuold,nyp,nx,inorder)
				! add guard cells for current in y direction
!				call paguard(cuold,kstrt,nvp,nx,kyp,ngds)	

				! Now find cuperp(t-dt)
				isign = -1
				call fft(cuold,cut,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
				call cuperp_jf(cut,nx,ny,kstrt)
				isign = 1
				call fft(cuold,cut,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
				! copy guard cells for current in y direction
				call pcguard(cuold,kstrt,nvp,kyp,inorder)
				! copy guard cells for current in x direction
				call cguard(cuold,nyp,nx,inorder)
				
				do icnt=1,size(cuold,2)
					do jcnt=1,size(cuold,3)
						ESPoynt_temp(1,icnt,jcnt,1) = cuold(1,icnt,jcnt,1) * phi_old(icnt,jcnt,1)
						ESPoynt_temp(2,icnt,jcnt,1) = cuold(2,icnt,jcnt,1) * phi_old(icnt,jcnt,1)
					enddo
				enddo
!						ESPoynt_temp(1,:,:,:) = cuold(1,:,:,:) * phi_old
!						ESPoynt_temp(2,:,:,:) = cuold(2,:,:,:) * phi_old
				!!!! Must mult by inv_num_par_cell when written out!
				! ESPoynt_temp now holds j_perp(t-dt)*phi(t-dt)*num_par_cell 



				ESPoynt_sum = ESPoynt_temp + ESPoynt_sum
				! copy guard cells for current in y direction
				call pcguard(ESPoynt_sum,kstrt,nvp,kyp,inorder)
				! copy guard cells for current in x direction
				call cguard(ESPoynt_sum,nyp,nx,inorder)


				
				isign = -1
				call fft(ESPoynt_temp,cut,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
				call ipdivf2(cut,phit,nx,ny,kstrt)				
				isign = 1
				call fft(div_ESPoynt,phit,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
				! copy guard cells for current in y direction
				call pcguard(div_ESPoynt,kstrt,nvp,kyp,inorder)
				! copy guard cells for current in x direction
				call cguard(div_ESPoynt,nyp,nx,inorder)
				div_ESPoynt_int = div_ESPoynt_int + div_ESPoynt















				cuold = cu !Store j(t-dt/2) for next time around
				
! transform potential to real space
				call pois(qt,sfieldt,ffc,temp_we,nx,ny,kstrt)
				isign = 1
				call fft(phi_old,sfieldt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
! copy to guard cells
				call pcguard(phi_old,kstrt,nvp,kyp,inorder)
				call cguard(phi_old,nyp,nx,inorder)		! phi_old = phi(t), which will be phi(t-dt) next time
																							! it's used.




				!To get vdote and poynt and u to be in same units, mult poynt and u by num_par_cell
				!not sure why need to divide by two
				if (nvdotE_int .ne. 0) then
					it = itime/nvdotE_int
					if (itime == nvdotE_int*it+1) then
						if (vdotEx_sumrec == 0) then
							! do x component
							vdotEx_sumrec = -1
							fname = './DIAG/vdotEx_int/vdotEx'
							vdotEx_sumfunit = get_funit(20)
							sfield = vdotE_sum(1,:,:,:)*dt/2.
							call writef(sfield,nxe,nypmx,vdotEx_sumfunit,-1,itime-1,itime*dt-dt,&
								&VDOTEX_INT,trim(fname),inorder)
							!do y component
							vdotEy_sumrec = -1
							fname = './DIAG/vdotEy_int/vdotEy'
							vdotEy_sumfunit = get_funit(20)
							sfield = vdotE_sum(2,:,:,:)*dt/2.
							call writef(sfield,nxe,nypmx,vdotEx_sumfunit,-1,itime-1,itime*dt-dt,&
								VDOTEY_INT,trim(fname),inorder)
						else
							sfield = vdotE_sum(1,:,:,:)*dt/2.
							fname = './DIAG/vdotEx_int/vdotEx'
							call writef(sfield,nxe,nypmx,vdotEx_sumfunit,-1,itime-1,itime*dt-dt,&
								&VDOTEX_INT,trim(fname),inorder)
							sfield = vdotE_sum(2,:,:,:)*dt/2.
							fname = './DIAG/vdotEy_int/vdotEy'
							call writef(sfield,nxe,nypmx,vdotEx_sumfunit,-1,itime-1,itime*dt-dt,&
								VDOTEY_INT,trim(fname),inorder)
						endif
					endif
				endif
				if (ntESPoynt_int > 0) then
					it = itime/ntESPoynt_int
					if (itime == ntESPoynt_int*it+1) then  !This is +1 because this diagnostic is 
						if (nESPoynt_sum_xrec == 0) then
							nESPoynt_sum_xrec = -1
							fname = './DIAG/ESPoynt_int_x/ESPoynt_int_x'
							ESPoynt_sum_xfunit = get_funit(20)
							sfield = ESPoynt_sum(1,:,:,:)*dt
							call writef(sfield,nxe,nypmx,ESPoynt_sum_xfunit,-1,itime-1,itime*dt-dt,&
								&ESPOYNT_INT_X,trim(fname),inorder)
						
							nESPoynt_sum_yrec = -1
							fname = './DIAG/ESPoynt_int_y/ESPoynt_int_y'
							ESPoynt_sum_yfunit = get_funit(20)
							sfield = ESPoynt_sum(2,:,:,:)*dt
							call writef(sfield,nxe,nypmx,ESPoynt_sum_yfunit,-1,itime-1,itime*dt-dt,&
								&ESPOYNT_INT_Y,trim(fname),inorder)
								
						else 
							fname = './DIAG/ESPoynt_int_x/ESPoynt_int_x'
							sfield = ESPoynt_sum(1,:,:,:)*dt
							call writef(sfield,nxe,nypmx,ESPoynt_sum_xfunit,-1,itime-1,itime*dt-dt,&
								&ESPOYNT_INT_X,trim(fname),inorder)
		
							fname = './DIAG/ESPoynt_int_y/ESPoynt_int_y'
							sfield = ESPoynt_sum(2,:,:,:)*dt
							call writef(sfield,nxe,nypmx,ESPoynt_sum_yfunit,-1,itime-1,itime*dt-dt,&
								&ESPOYNT_INT_Y,trim(fname),inorder)
						endif
					endif
				endif
				if (ntfield_ene_int > 0) then
					it = itime/ntfield_ene_int
					if (itime==ntfield_ene_int*it+1) then
						if (ntfield_ene_int_xrec==0) then
							ntfield_ene_int_xrec = -1
							ntfield_ene_int_yrec = -1
							fname = './DIAG/u_x_int/u_x_int'
							sfield = grad_phi_dt(1,:,:,:)*num_par_cell
							ntfield_ene_int_funitx = get_funit(20)
							call writef(sfield,nxe,nypmx,ntfield_ene_int_funitx,-1,itime-1,itime*dt-dt,&
								&EX_ENE_INT,trim(fname),inorder)
							fname = './DIAG/u_y_int/u_y_int'
							sfield = grad_phi_dt(2,:,:,:)*num_par_cell
							ntfield_ene_int_funity = get_funit(20)
							call writef(sfield,nxe,nypmx,ntfield_ene_int_funity,-1,itime-1,itime*dt-dt,&
								&EY_ENE_INT,trim(fname),inorder)
						else
							sfield = grad_phi_dt(1,:,:,:)*num_par_cell
							fname = './DIAG/u_x_int/u_x_int'
							call writef(sfield,nxe,nypmx,ntfield_ene_int_funitx,-1,itime-1,itime*dt-dt,&
								&EX_ENE_INT,trim(fname),inorder)
							fname = './DIAG/u_y_int/u_y_int'
							sfield = grad_phi_dt(2,:,:,:)*num_par_cell
							call writef(sfield,nxe,nypmx,ntfield_ene_int_funity,-1,itime-1,itime*dt-dt,&
								&EY_ENE_INT,trim(fname),inorder)
						endif
					endif
				endif
				if (nt_div_ESPoynt_int > 0) then
					it = itime/nt_div_ESPoynt_int
					if (itime==nt_div_ESPoynt_int*it+1) then
						if (div_ESPoynt_int_rec==0) then
							div_ESPoynt_int_rec = -1
							fname = './DIAG/div_ESPoynt_int/div_ESPoynt_int'
							sfield = div_ESPoynt_int(:,:,:)*dt
							div_ESPoynt_int_funit = get_funit(20)
							call writef(sfield,nxe,nypmx,div_ESPoynt_int_funit,-1,itime-1,itime*dt-dt,&
								&DIV_ESPOYNTINT,trim(fname),inorder)
						else
							fname = './DIAG/div_ESPoynt_int/div_ESPoynt_int'
							sfield = div_ESPoynt_int(:,:,:)*dt
							call writef(sfield,nxe,nypmx,div_ESPoynt_int_funit,-1,itime-1,itime*dt-dt,&
								&DIV_ESPOYNTINT,trim(fname),inorder)
						endif
					endif
				endif
			endif
						
			ES_ene = 0.
			vdote_tot = 0.
			poynt = 0.
			div_poynt = 0.
			top=0.
			bot=0.
			left=0.
			right=0.
			do icnt=2,31
				do jcnt=2,ny+1
!			do icnt=16,20
!				do jcnt=16,20
!			do icnt=2,nx
!				do jcnt=2,ny+1
					
					top = top + 0.5*(grad_phi(1,icnt,jcnt,1)**2+grad_phi(2,icnt,jcnt,1)**2)

					ES_ene = ES_ene + 0.5*( grad_phi(1,icnt,jcnt,1)**2 + grad_phi(2,icnt,jcnt,1)**2 - &
																	(grad_phi_dt2(1,icnt,jcnt,1)**2 + grad_phi_dt2(2,icnt,jcnt,1)**2) )/2./dt
!					ES_ene = ES_ene + grad_phi_dt(1,icnt,jcnt,1) + grad_phi_dt(2,icnt,jcnt,1)

					vdote_tot = vdote_tot + vdotE_inst(1,icnt,jcnt,1) + vdotE_inst(2,icnt,jcnt,1)
!					vdote_tot = vdote_tot + vdotE_sum(1,icnt,jcnt,1) + vdotE_sum(2,icnt,jcnt,1)
					
					div_poynt = div_poynt + div_ESPoynt(icnt,jcnt,1)		!Using Viktor's divergence routine and putting it into phias temp
!					div_poynt = div_poynt + div_ESPoynt_int(icnt,jcnt,1)		!Using Viktor's divergence routine and putting it into phias temp
!					poynt = poynt + ESPoynt_sum(1,icnt,jcnt,1) + ESPoynt_sum(2,icnt,jcnt,1)
				enddo
			enddo
!			do icnt=2,nx+1
!				top = top + ESPoynt_sum(2,icnt,ny+1,1)	!top
!				bot = bot - ESPoynt_sum(2,icnt,1,1)		!bot
!			enddo
!			do icnt=2,ny+1
!				right = right + ESPoynt_sum(1,nx+1,icnt,1)	!right
!				left = left - ESPoynt_sum(1,1,icnt,1)		!left
!			enddo
!			do icnt=2,21
!				top = top + ESPoynt_sum(2,icnt,21,1)	!top
!				bot = bot - ESPoynt_sum(2,icnt,2,1)		!bot
!				right = right - ESPoynt_sum(1,21,icnt,1)	!right
!				left = left + ESPoynt_sum(1,2,icnt,1)		!left
!			enddo
!			poynt = (top+bot+left+right)*dt
			
!			allthree = ES_ene*num_par_cell+vdote_tot/2.*dt+div_poynt
!			print*,itime,ES_ene*num_par_cell,vdote_tot/2.*dt,div_poynt,allthree

!			print*,itime,ES_ene*num_par_cell,vdote_tot,div_poynt,ES_ene*num_par_cell+vdote_tot+div_poynt,top*num_par_cell,wtot(1)



!			print*,poynt,div_poynt*dt

!			print*,itime,ES_ene*num_par_cell,vdote_tot/2.,wtot
!			print*,itime, top,right,bot,left,poynt
			
!			kinene = 0.
!			do icnt=1, npp(1)
!				if ((part(1,icnt,1) >1.) .and. (part(1,icnt,1) < 19.) ) then
!				if ((part(2,icnt,1) >1.) .and. (part(2,icnt,1) < 19.) ) then
!					kinene = kinene + 0.5*(part(3,icnt,1)**2 + part(4,icnt,1)**2)
!				endif
!				endif
!			enddo
!			print*,kinene-kinene_old,vdote_tot*inv_num_par_cell
!			print*,kinene-kinene_old,vdote_tot*dt
!			kinene_old2 = kinene_old
!			if (itime == 1) then
!				wtotold = wtot(4)
!				wkeold = wtot(2)
!				kinene_old = kinene
!			endif
			
!******************************************************************************************************
