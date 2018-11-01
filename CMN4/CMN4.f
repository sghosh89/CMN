! allocated carbon C = x, mutualist M = y, non-mutualist N = z, construction carbon = v
          
             implicit real*8(a-h,o-z)
                common/par1/ps  !P availability in soil
                common/par2/akc !half-saturation constant of allocated C
                common/par3/u    !rate of P return by mutualist
                common/par4/af   !fidelity of plant allocation to mutualist symbiont
                common/par5/bmax !maximum growth rate of symbionts
                common/par6/as  !cost of mutualism
                common/par7/d !death rate of symbionts
                common/par8/akm !half-saturation constant for mutualist
                common/par9/akn !half-saturation constant for non-mutualist
                common/par10/phi  !construction C initial value
                common/par11/g  !construction C decay rate
                
                double precision ps,akc,u,af
                double precision bmax,as,d,akm
                double precision akn,phi,g
             
             open(10,file='CMN4.in',status='unknown')
             open(20,file='CMN4_tCM.dat',status='unknown')
             open(30,file='CMN4_tCN.dat',status='unknown')
             open(40,file='CMN4_tMN.dat',status='unknown')
             open(50,file='CMN4_tCaCc.dat',status='unknown')
             open(60,file='CMN4_phiCrSr.dat',status='unknown')
                                
             read(10,*)ps,akc,u,af
             read(10,*)bmax,as,d,akm
             !read(10,*)akn,phi,g
             read(10,*)akn,g
             read(10,*)hh,kk
             read(10,*)x0,y0,z0,v0,t0
            
             as_max=(1.0-(d/bmax))                              !condtion for coexistence 
             print *, "as_max=",as_max
             print *,"as=",as
              Cam=(akm*d)/((bmax*(1.0-as))-d)
              Can=(akn*d)/((bmax-d)*((1.0-af)**2.0))
              print *,"Can*/Cam*=",Can/Cam
              print *,"1/(1-f)=",1/(1.0-af)
          
          phi0=3.0
          
          do ip=0,294
             phi=phi0+(0.5*ip)
             print *,"----------phi=",phi,"--------------"
     
          do i=1,kk
            
              ak1f1=hh*f1(x0,y0,z0,t0)
              ak1f2=hh*f2(x0,y0,v0,t0)
              ak1f3=hh*f3(x0,z0,v0,t0)
              ak1f4=hh*f4(y0,z0,v0,t0)
                   
                    temp1f1=x0+0.5*ak1f1
                    temp1f2=y0+0.5*ak1f2
                    temp1f3=z0+0.5*ak1f3
                    temp1f4=v0+0.5*ak1f4

             ak2f1=hh*f1(temp1f1,temp1f2,temp1f3,t0+0.5*hh)
             ak2f2=hh*f2(temp1f1,temp1f2,temp1f4,t0+0.5*hh)
             ak2f3=hh*f3(temp1f1,temp1f3,temp1f4,t0+0.5*hh)
             ak2f4=hh*f4(temp1f2,temp1f3,temp1f4,t0+0.5*hh)
            
                    temp2f1=x0+0.5*ak1f1
                    temp2f2=y0+0.5*ak1f2
                    temp2f3=z0+0.5*ak1f3
                    temp2f4=v0+0.5*ak1f4

             ak3f1=hh*f1(temp2f1,temp2f2,temp2f3,t0+0.5*hh)
             ak3f2=hh*f2(temp2f1,temp2f2,temp2f4,t0+0.5*hh)
             ak3f3=hh*f3(temp2f1,temp2f3,temp2f4,t0+0.5*hh)
             ak3f4=hh*f4(temp2f2,temp2f3,temp2f4,t0+0.5*hh)
             
                    temp3f1=x0+ak3f1
                    temp3f2=y0+ak3f2
                    temp3f3=z0+ak3f3
                    temp3f4=v0+ak3f4

             ak4f1=hh*f1(temp3f1,temp3f2,temp3f3,t0+hh)
             ak4f2=hh*f2(temp3f1,temp3f2,temp3f4,t0+hh)
             ak4f3=hh*f3(temp3f1,temp3f3,temp3f4,t0+hh)
             ak4f4=hh*f4(temp3f2,temp3f3,temp3f4,t0+hh)

      
                     x=x0+(ak1f1+2.0*(ak2f1+ak3f1)+ak4f1)/6.0
                     y=y0+(ak1f2+2.0*(ak2f2+ak3f2)+ak4f2)/6.0
                     z=z0+(ak1f3+2.0*(ak2f3+ak3f3)+ak4f3)/6.0
                     v=v0+(ak1f4+2.0*(ak2f4+ak3f4)+ak4f4)/6.0

             write(20,*)t0,x0,y0
             write(30,*)t0,x0,z0
             write(40,*)t0,y0,z0
             write(50,*)t0,x0,v0
 
        
               t0=t0+hh
               x0=x
               y0=y
               z0=z
               v0=v

             end do

             print *,"Ca=",x
             print *,"Cc=",v
             print *,"M=",y
             print *,"N=",z

           Cr=v/x
           Sr=z/y
           write(60,*)phi,Cr,Sr
            
          end do

!             stop
             end

!----------------------------------------------------------------------------------------------------------

             function f1(x,y,z,t)
        
             implicit real*8(a-h,o-z)
                common/par1/ps  !P availability in soil
                common/par2/akc !half-saturation constant of allocated C
                common/par3/u    !rate of P return by mutualist
                common/par4/af   !fidelity of plant allocation to mutualist symbiont
                common/par5/bmax !maximum growth rate of symbionts
                common/par6/as  !cost of mutualism
                common/par7/d !death rate of symbionts
                common/par8/akm !half-saturation constant for mutualist
                common/par9/akn !half-saturation constant for non-mutualist
                common/par10/phi  !construction C
                common/par11/g  !construction C decay rate
          f1=(x*(y/(akc+y))*u)*((y/(y+z))/(1.0-akf+(akf*(y/(y+z)))))  
          f1=1.0-ps-f1

                    return
                    end
!-------------------------------------------------------------------------------------------------------------------

              function f2(x,y,v,t)

             implicit real*8(a-h,o-z)
                common/par1/ps  !P availability in soil
                common/par2/akc !half-saturation constant of allocated C
                common/par3/u    !rate of P return by mutualist
                common/par4/af   !fidelity of plant allocation to mutualist symbiont
                common/par5/bmax !maximum growth rate of symbionts
                common/par6/as  !cost of mutualism
                common/par7/d !death rate of symbionts
                common/par8/akm !half-saturation constant for mutualist
                common/par9/akn !half-saturation constant for non-mutualist
                common/par10/phi  !construction C
                common/par11/g  !construction C decay rate
               
        f2=((bmax*(1.0-as)*(x+v)/(akm+x+v))-d)*y

                return
                end
!------------------------------------------------------------------------------------------------------------------------

            function f3(x,z,v,t)

             implicit real*8(a-h,o-z)
                common/par1/ps  !P availability in soil
                common/par2/akc !half-saturation constant of allocated C
                common/par3/u    !rate of P return by mutualist
                common/par4/af   !fidelity of plant allocation to mutualist symbiont
                common/par5/bmax !maximum growth rate of symbionts
                common/par6/as  !cost of mutualism
                common/par7/d !death rate of symbionts
                common/par8/akm !half-saturation constant for mutualist
                common/par9/akn !half-saturation constant for non-mutualist
                common/par10/phi  !construction C
                common/par11/g  !construction C decay rate
               
           f3=bmax*(1.0-af)*((x*(1.0-af))+v)
           f3=(f3/(akn+((1.0-af)*((x*(1.0-af))+v))))-d
           f3=(f3*z)
            
               return
               end

!------------------------------------------------------------------------------------------------------------------------

            function f4(y,z,v,t)

             implicit real*8(a-h,o-z)
                common/par1/ps  !P availability in soil
                common/par2/akc !half-saturation constant of allocated C
                common/par3/u    !rate of P return by mutualist
                common/par4/af   !fidelity of plant allocation to mutualist symbiont
                common/par5/bmax !maximum growth rate of symbionts
                common/par6/as  !cost of mutualism
                common/par7/d !death rate of symbionts
                common/par8/akm !half-saturation constant for mutualist
                common/par9/akn !half-saturation constant for non-mutualist
                common/par10/phi  !construction C
                common/par11/g  !construction C decay rate
               
           f4=phi-(g*v*(y+z))
               return
               end










            
