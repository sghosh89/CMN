! 4 state variables: 
! 	Carbon allocation rate by the plant = x, 
! 	Density of mutualist symbiont = y, 
!	Density of non-mutualist symbiont = z, 
!	Length of new roots available for colonization by symbionts = v

! Consider for these default parameters in the input file:
!	Half saturation constant for the symbionts aKM = aKN = 10, 
!	Half saturation constant for allocated carbon aKA =5, 
!	Energy allocation rate by the plant per colonized root eM = eN = 0.5, 
!	Phosphorus uptake per unit of preferentially allocated carbon received by the mutualists u=0.4, 
!	Maximum growth rate of symbionts per capita bmax=0.8, 
!	Death rate of symbionts per capita d=0.5
          
             implicit real*8(a-h,o-z)
               
                common/par1/ps  !concentration of soil phosphorous available for plant
                common/par2/u      !phosphorous uptake per unit of A ??unit
                common/par3/aKA    !half saturation constant for A
                common/par4/af     ! fidelity of plant allocation to mutualist
          	common/par5/phi     ! base root growth rate
          	common/par6/aM     ! colonization rate of new roots by M
          	common/par7/aN     ! colonization rate of new roots by N
          	common/par8/eM     ! plant's efficiency to give resources to mutualist
          	common/par9/eN     !and non-mutualist eM should be equal to eN
          	common/par10/bmax   ! max growth rate for symbionts
          	common/par11/s      ! cost of mutualism
          	common/par12/aKM     ! half saturation constant for M
          	common/par13/d      ! death rate for symbionts
          	common/par14/aKN     ! half saturation constant for N
                
                double precision ps,u,aKA,af,ps0
                double precision phi,aM,aN,eM
                double precision eN,bmax,s,aKM
                double precision d,aKN
             
             open(10,file='ARMN.in',status='unknown')
             open(40,file='ARMN_tAR.dat',status='unknown')
             open(50,file='ARMN_tMN.dat',status='unknown')
 
                                             
             read(10,*)ps,u,aKA,af
             read(10,*)phi,aM,aN
             read(10,*)eM,eN,bmax,d
             read(10,*)s,aKM,aKN
             read(10,*)hh,kk
             read(10,*)x0,y0,z0,v0,t0

             af_min=1.0-((aKN/aKM)*(((bmax*(1.0-s))-d)/(bmax-d)))
 
             print *,"----------given f=",af,"--------------"
             print *,"----------f_min=",af_min,"--------------"
            
    

     
            do i=1,kk
            
              ak1f1=hh*f1(x0,y0,z0,t0)
              ak1f2=hh*f2(x0,y0,z0,v0,t0)
              ak1f3=hh*f3(x0,y0,z0,v0,t0)
              ak1f4=hh*f4(y0,z0,v0,t0)
                   
                    temp1f1=x0+0.5*ak1f1
                    temp1f2=y0+0.5*ak1f2
                    temp1f3=z0+0.5*ak1f3
                    temp1f4=v0+0.5*ak1f4

             ak2f1=hh*f1(temp1f1,temp1f2,temp1f3,t0+0.5*hh)
             ak2f2=hh*f2(temp1f1,temp1f2,temp1f3,temp1f4,t0+0.5*hh)
             ak2f3=hh*f3(temp1f1,temp1f2,temp1f3,temp1f4,t0+0.5*hh)
             ak2f4=hh*f4(temp1f2,temp1f3,temp1f4,t0+0.5*hh)
            
                    temp2f1=x0+0.5*ak1f1
                    temp2f2=y0+0.5*ak1f2
                    temp2f3=z0+0.5*ak1f3
                    temp2f4=v0+0.5*ak1f4

             ak3f1=hh*f1(temp2f1,temp2f2,temp2f3,t0+0.5*hh)
             ak3f2=hh*f2(temp2f1,temp2f2,temp2f3,temp2f4,t0+0.5*hh)
             ak3f3=hh*f3(temp2f1,temp2f2,temp2f3,temp2f4,t0+0.5*hh)
             ak3f4=hh*f4(temp2f2,temp2f3,temp2f4,t0+0.5*hh)
             
                    temp3f1=x0+ak3f1
                    temp3f2=y0+ak3f2
                    temp3f3=z0+ak3f3
                    temp3f4=v0+ak3f4

             ak4f1=hh*f1(temp3f1,temp3f2,temp3f3,t0+hh)
             ak4f2=hh*f2(temp3f1,temp3f2,temp3f3,temp3f4,t0+hh)
             ak4f3=hh*f3(temp3f1,temp3f2,temp3f3,temp3f4,t0+hh)
             ak4f4=hh*f4(temp3f2,temp3f3,temp3f4,t0+hh)

      
                     x=x0+(ak1f1+2.0*(ak2f1+ak3f1)+ak4f1)/6.0
                     y=y0+(ak1f2+2.0*(ak2f2+ak3f2)+ak4f2)/6.0
                     z=z0+(ak1f3+2.0*(ak2f3+ak3f3)+ak4f3)/6.0
                     v=v0+(ak1f4+2.0*(ak2f4+ak3f4)+ak4f4)/6.0

             write(40,*)t0,x0,v0
             write(50,*)t0,y0,z0
 
        
               t0=t0+hh
               x0=x
               y0=y
               z0=z
               v0=v

             end do

             print *,"Caeq=",x
             print *,"Cceq=",v
             print *,"Meq=",y
             print *,"Neq=",z

          !evaluate at equilibrium
          f1eqm=f1(x0,y0,z0,t0)
          f2eqm=f2(x0,y0,z0,v0,t0)
          f3eqm=f3(x0,y0,z0,v0,t0)
          f4eqm=f4(y0,z0,v0,t0)

         print *,"derv_Ca_ateqm=",f1eqm
         print *,"derv_Cc_ateqm=",f4eqm
         print *,"derv_M_ateqm=",f2eqm
         print *,"derv_N_ateqm=",f3eqm

             stop
             end

!----------------------------------------------------------------------------------------------------------

             function f1(x,y,z,t)
        
             implicit real*8(a-h,o-z)
                 common/par1/ps  !concentration of soil phosphorous available for plant
                common/par2/u      !phosphorous uptake per unit of A ??unit
                common/par3/aKA    !half saturation constant for A
                common/par4/af     ! fidelity of plant allocation to mutualist
          	common/par5/phi     ! base root growth rate
          	common/par6/aM     ! colonization rate of new roots by M
          	common/par7/aN     ! colonization rate of new roots by N
          	common/par8/eM     ! plant's efficiency to give resources to mutualist
          	common/par9/eN     !and non-mutualist eM should be equal to eN
          	common/par10/bmax   ! max growth rate for symbionts
          	common/par11/s      ! cost of mutualism
          	common/par12/aKM     ! half saturation constant for M
          	common/par13/d      ! death rate for symbionts
          	common/par14/aKN     ! half saturation constant for N

          	pM= y/(y+z) !proportion of mutualist
         	FMN= u*(y/(aKA+y))*(pM/(1-af+(af*pM))) !saturating function  
           
  	        f1=1.0-ps-(FMN*x)

                return
                end
!-------------------------------------------------------------------------------------------------------------------

              function f2(x,y,z,v,t)

             implicit real*8(a-h,o-z)
                common/par1/ps  !concentration of soil phosphorous available for plant
                common/par2/u      !phosphorous uptake per unit of A ??unit
                common/par3/aKA    !half saturation constant for A
                common/par4/af     ! fidelity of plant allocation to mutualist
          	common/par5/phi     ! base root growth rate
          	common/par6/aM     ! colonization rate of new roots by M
          	common/par7/aN     ! colonization rate of new roots by N
          	common/par8/eM     ! plant's efficiency to give resources to mutualist
          	common/par9/eN     !and non-mutualist eM should be equal to eN
          	common/par10/bmax   ! max growth rate for symbionts
          	common/par11/s      ! cost of mutualism
          	common/par12/aKM     ! half saturation constant for M
          	common/par13/d      ! death rate for symbionts
          	common/par14/aKN     ! half saturation constant for N
               
                pM= y/(y+z) !proportion of mutualist
                denom=(y+z)*(1.0-af+(af*pM))
                CM= (eM*aM*v) + (x/denom) !carbon consumed by Mutualist    
               
                f2=(((bmax*(1.0-s)*CM)/(aKM+CM))-d)*y

                return
                end
!------------------------------------------------------------------------------------------------------------------------

            function f3(x,y,z,v,t)

             implicit real*8(a-h,o-z)
                 common/par1/ps  !concentration of soil phosphorous available for plant
                common/par2/u      !phosphorous uptake per unit of A ??unit
                common/par3/aKA    !half saturation constant for A
                common/par4/af     ! fidelity of plant allocation to mutualist
          	common/par5/phi     ! base root growth rate
          	common/par6/aM     ! colonization rate of new roots by M
          	common/par7/aN     ! colonization rate of new roots by N
          	common/par8/eM     ! plant's efficiency to give resources to mutualist
          	common/par9/eN     !and non-mutualist eM should be equal to eN
          	common/par10/bmax   ! max growth rate for symbionts
          	common/par11/s      ! cost of mutualism
          	common/par12/aKM     ! half saturation constant for M
          	common/par13/d      ! death rate for symbionts
          	common/par14/aKN     ! half saturation constant for N
                
               pM= y/(y+z) !proportion of mutualist
               denom=(y+z)*(1.0-af+(af*pM))
               CN= (eN*aN*v) + ((x*(1.0-af))/denom) !carbon consumed by Non-mutualist
               f3=(((bmax*CN)/(aKN+CN))-d)*z
            
               return
               end

!------------------------------------------------------------------------------------------------------------------------

            function f4(y,z,v,t)

             implicit real*8(a-h,o-z)
                 common/par1/ps  !concentration of soil phosphorous available for plant
                common/par2/u      !phosphorous uptake per unit of A ??unit
                common/par3/aKA    !half saturation constant for A
                common/par4/af     ! fidelity of plant allocation to mutualist
          	common/par5/phi     ! base root growth rate
          	common/par6/aM     ! colonization rate of new roots by M
          	common/par7/aN     ! colonization rate of new roots by N
          	common/par8/eM     ! plant's efficiency to give resources to mutualist
          	common/par9/eN     !and non-mutualist eM should be equal to eN
          	common/par10/bmax   ! max growth rate for symbionts
          	common/par11/s      ! cost of mutualism
          	common/par12/aKM     ! half saturation constant for M
          	common/par13/d      ! death rate for symbionts
          	common/par14/aKN     ! half saturation constant for N
                
               
               f4=phi-(aM*v*y)-(aN*v*z)
               return
               end










            
