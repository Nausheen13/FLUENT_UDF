#include "udf.h"
#include "dpm.h"
#include "stdio.h"
#include "sg_mphase.h"
#include "sg.h"
#include "sg_mem.h"
#include "mem.h"
#include "flow.h"

#define Fluid_mass 0.000444



DEFINE_SOURCE(zMomentum, c,t,dS,eqn)
{double A,B,fraction,ConstantV,ConstantI;
 double C1=931.684e+03, C2=584.505e+08;
 real mass, Uz, source;
 
 
 
 A= C_DPMS_CONCENTRATION(c,t)*C_VOLUME(c,t)*46110;
 B= A*(CURRENT_TIME);
 
 
 mass= C_R(c,t)*C_VOLUME(c,t);
 fraction= B/(0.000444+B);


 
 ConstantI= C1+(C1*fraction);
 ConstantV= C2+(C2*fraction);
 
 
 
 Uz= C_W(c,t);
 source= -((C_MU_L(c,t)*ConstantV*Uz)+(0.5*ConstantI*C_R(c,t)*fabs(Uz)*Uz));
 dS[eqn]= -((C_MU_L(c,t)*ConstantV)+(2*ConstantI*C_R(c,t)*fabs(Uz)));
 
 return source;
}
 
DEFINE_SOURCE(xMomentum, c,t,dS,eqn)
{double A,B,fraction,ConstantV,ConstantI;
 double C1=931.684e+03, C2=584.505e+08;
 real mass, Ux, source;
 
 A= C_DPMS_CONCENTRATION(c,t)*C_VOLUME(c,t)*46110;
 B= A*(CURRENT_TIME);
 
 mass= C_R(c,t)*C_VOLUME(c,t);
 fraction= B/(0.000444+B);
 
 ConstantI= C1+(C1*fraction);
 ConstantV= C2+(C2*fraction);
 
 Ux= C_U(c,t);
 source= -((C_MU_L(c,t)*ConstantV*Ux)+(0.5*ConstantI*C_R(c,t)*fabs(Ux)*Ux));
 dS[eqn]= -((C_MU_L(c,t)*ConstantV)+(2*ConstantI*C_R(c,t)*fabs(Ux)));
 
 return source;
}
 
 DEFINE_SOURCE(yMomentum, c,t,dS,eqn)
{double A,B,fraction,ConstantV,ConstantI;
 double C1=931.684e+03, C2=584.505e+08;
 real mass, Uy, source;
 
 A= C_DPMS_CONCENTRATION(c,t)*C_VOLUME(c,t)*46110;
 B= A*(CURRENT_TIME);
 
 mass= C_R(c,t)*C_VOLUME(c,t);
 fraction= B/(0.000444+B);
 
 ConstantI= C1+(C1*fraction);
 ConstantV= C2+(C2*fraction);
 
 C_UDMI(c,t,2)= ConstantI;
 C_UDMI(c,t,3)= ConstantV;
 
 Uy= C_V(c,t);
 source= -((C_MU_L(c,t)*ConstantV*Uy)+(0.5*ConstantI*C_R(c,t)*fabs(Uy)*Uy));
 dS[eqn]= -((C_MU_L(c,t)*ConstantV)+(2*ConstantI*C_R(c,t)*fabs(Uy)));
 

 
 return source;
}

Define_DELTAT(mydeltat,domain)
{real T,time_step;
 T= CURRENT_TIME;
 
 if(T>300)
 {
  time_step= (CURRENT_TIMESTEP*10);
  }
  
  else
  {
   time_step= CURRENT_TIMESTEP;
  }
  
  return time_step;
 }
