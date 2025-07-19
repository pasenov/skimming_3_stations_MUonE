#ifndef Station_hits_H
#define Station_hits_H

#include <array>

#include <iostream>

using namespace std;


struct Station_hits {
  int fired_xy_modules = 0;
  int multifired_xy_modules = 0;
  int threefired_xy_modules = 0;
  int fourfired_xy_modules = 0;
  int fivefired_xy_modules = 0;
  int nhits_xy = 0;
  int fired_uv_modules = 0;
  int multifired_uv_modules = 0;
  int threefired_uv_modules = 0;
  int fourfired_uv_modules = 0;
  int fivefired_uv_modules = 0;
  int nhits_uv = 0;
  //additional condition: more than one hit in each of last modules 
  int fired_last_modules = 0;
	int nhits_last = 0;
	int multifired_last_modules = 0;
	int threefired_last_modules = 0;
	int fourfired_last_modules = 0;
  int fivefired_last_modules = 0;

  int fired_first_modules = 0;
	int nhits_first = 0;
	int multifired_first_modules = 0;
	int threefired_first_modules = 0;
	int fourfired_first_modules = 0;
  int fivefired_first_modules = 0;

};

Station_hits fill_station_hits(const std::array<int,18> & stubs, int ilay0) {
  
  Station_hits S;
  
  for (int jm=0; jm<6; ++jm) {
    
    int ns = stubs[ilay0+jm];
      
    if (ns>0)
      {
	if (jm==2 || jm==3)
	  {
	    S.fired_uv_modules++;
	    S.nhits_uv = S.nhits_uv + ns;
	    if (ns>1) S.multifired_uv_modules++;
	    if (ns>2) S.threefired_uv_modules++;
	    if (ns>3) S.fourfired_uv_modules++;
	    if (ns>4) S.fivefired_uv_modules++;
	  }
	else
	  {
	    S.fired_xy_modules++;
	    S.nhits_xy = S.nhits_xy + ns;
	    if (ns>1) S.multifired_xy_modules++;
	    if (ns>2) S.threefired_xy_modules++;
	    if (ns>3) S.fourfired_xy_modules++;
	    if (ns>4) S.fivefired_xy_modules++;
      //additional condition
      if (jm==4 || jm==5)
	    {
	    S.fired_last_modules++;
	    S.nhits_last = S.nhits_last + ns;
	    if (ns>1) S.multifired_last_modules++;
	    if (ns>2) S.threefired_last_modules++;
	    if (ns>3) S.fourfired_last_modules++;
	    if (ns>4) S.fivefired_last_modules++;
	    }
      else if (jm==0 || jm==1)
	    {
	    S.fired_first_modules++;
	    S.nhits_first = S.nhits_first + ns;
	    if (ns>1) S.multifired_first_modules++;
	    if (ns>2) S.threefired_first_modules++;
	    if (ns>3) S.fourfired_first_modules++;
	    if (ns>4) S.fivefired_first_modules++;
	    }

    
	  }
      }
  }

  // if (S.multifired_xy_modules > 0 && S.multifired_uv_modules > 0) {
  
  //   cout<<".................................................."<<endl
  //   <<" fired_xy_modules = "<< S.fired_xy_modules
  //   <<" multifired_xy_modules = "<< S.multifired_xy_modules
  //   <<" threefired_xy_modules = "<< S.threefired_xy_modules
  //   <<" fourfired_xy_modules = "<< S.fourfired_xy_modules
  //   <<" fivefired_xy_modules = "<< S.fivefired_xy_modules
  //   <<" nhits_xy = "<< S.nhits_xy
  //   <<" fired_uv_modules = " << S.fired_uv_modules 
  //   <<" multifired_uv_modules = " << S.multifired_uv_modules 
  //   <<" threefired_uv_modules = "<< S.threefired_uv_modules
  //   <<" fourfired_uv_modules = "<< S.fourfired_uv_modules
  //   <<" fivefired_uv_modules = "<< S.fivefired_uv_modules
  //   <<" nhits_uv = " << S.nhits_uv
  //   <<".................................................."<<endl;

  // }
  
  
  return S;
};

#endif
