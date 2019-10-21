! Allocated carbon = x, Mutualist = y, Non-mutualist N = z, uncolonized Root = v
          
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
             
            open(10,file='ARMN_fps.in',status='unknown')

            open(88,file='ARMN_fpsAR.dat',status='unknown') ! see variation of A,R at eqm against fidelity and soil P
            open(90,file='ARMN_fpsMN.dat',status='unknown') ! see variation of M,N at eqm against fidelity and soil P
            open(92,file='ARMN_fpsPM.dat',status='unknown') ! see variation of PM(proportion of mutualist) at eqm against fidelity and soil P

                                             
             read(10,*)u,aKA             
             read(10,*)phi,aM,aN
             read(10,*)eM,eN,bmax,d
             read(10,*)s,aKM,aKN
             read(10,*)hh,kk

             af_min=1.0-((aKN/aKM)*(((bmax*(1.0-s))-d)/(bmax-d)))
             ps0=0.0
      
          
             af_incr=(1.0-af_min)/20.0 
             ps_incr=(1.0-ps0)/20.0

          do iaf=0,20
             af=af_min+(af_incr*iaf)
             
           do ip=0,20
             ps=ps0+(ps_incr*ip)
             
             !print *,"-----f=",af,"-----ps =",ps,"----"

               x0= 0.5
               y0= 0.1
               z0= 0.1
               v0= 0.5
               t0= 0.0      

     
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
        
               t0=t0+hh
               x0=x
               y0=y
               z0=z
               v0=v

             end do

             !print *,"Caeq=",x
             !print *,"Cceq=",v
             !print *,"Meq=",y
             !print *,"Neq=",z

          prop_y=y/(y+z) !proportion of mutualist
          
           write(88,*)af,ps,x,v
           write(90,*)af,ps,y,z           
           write(92,*)af,ps,prop_y

            
         end do
        end do

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










            
