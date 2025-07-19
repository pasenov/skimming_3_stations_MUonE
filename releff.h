#ifndef releff_H
#define releff_H

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "TRint.h"

std::pair<double,double> binomial_eff(Long64_t num, Long64_t den) {

  if (den !=0) {
    double eff = double(num)/double(den);
    double err = sqrt(eff*(1-eff)/double(den));
    return std::make_pair(eff, err);
  }
  else return {0,0};

};

void read_counts(std::ifstream & sum,
 Long64_t & den_all,
 Long64_t & den_presel, Long64_t & num_trackable, Long64_t & num_golden, Long64_t & num_zero,
 Long64_t & den_2presel, Long64_t & num_2trackable, Long64_t & num_2golden,
 Long64_t & den_3presel, Long64_t & num_3trackable, Long64_t & num_3golden,
 Long64_t & den_4presel, Long64_t & num_4trackable, Long64_t & num_4golden)
{
  std::string s, c;
  std::istringstream ss;
  //---
  getline(sum, s);
  //---
  getline(sum, s); ss = std::istringstream(s); ss>>c>>den_all;
  getline(sum, s); ss = std::istringstream(s); ss>>c>>den_presel;
  getline(sum, s); ss = std::istringstream(s); ss>>c>>num_trackable;
  getline(sum, s); ss = std::istringstream(s); ss>>c>>num_golden;
  //---  
  getline(sum, s);
  getline(sum, s);
  getline(sum, s);
  //---
  getline(sum, s); ss = std::istringstream(s); ss>>c>>num_zero;
  //---  
  getline(sum, s);
  //---  
  getline(sum, s); ss = std::istringstream(s); ss>>c>>den_2presel;
  getline(sum, s); ss = std::istringstream(s); ss>>c>>num_2trackable;
  getline(sum, s); ss = std::istringstream(s); ss>>c>>num_2golden;
  //---  
  getline(sum, s);
  getline(sum, s);
  //---  
  getline(sum, s); ss = std::istringstream(s); ss>>c>>den_3presel;
  getline(sum, s); ss = std::istringstream(s); ss>>c>>num_3trackable;
  getline(sum, s); ss = std::istringstream(s); ss>>c>>num_3golden;
  //---  
  getline(sum, s);
  getline(sum, s);
  //---  
  getline(sum, s); ss = std::istringstream(s); ss>>c>>den_4presel;
  getline(sum, s); ss = std::istringstream(s); ss>>c>>num_4trackable;
  getline(sum, s); ss = std::istringstream(s); ss>>c>>num_4golden;
  //---  
  getline(sum, s);
  getline(sum, s);
}

void add_efficiency(char inum, char iden, std::ifstream & sumold, std::ifstream & i, std::ofstream & sum) {

  Long64_t Sum_den_all,
    Sum_den_presel, Sum_num_trackable, Sum_num_golden,
    Sum_num_zero,
    Sum_den_2presel, Sum_num_2trackable, Sum_num_2golden,
    Sum_den_3presel, Sum_num_3trackable, Sum_num_3golden,
    Sum_den_4presel, Sum_num_4trackable, Sum_num_4golden;

  read_counts(sumold,
	      Sum_den_all,
	      Sum_den_presel, Sum_num_trackable, Sum_num_golden,
	      Sum_num_zero,
	      Sum_den_2presel, Sum_num_2trackable, Sum_num_2golden,
	      Sum_den_3presel, Sum_num_3trackable, Sum_num_3golden,
	      Sum_den_4presel, Sum_num_4trackable, Sum_num_4golden);
  
  std::cout<<"RELATIVE EFFICIENCY Station S_"<<inum
  	   <<" (in events with Golden Passing Mu in S_"<<iden<<")"<<std::endl;
  std::cout<<"sumold (1mu): "<<Sum_den_all<<", "<<Sum_den_presel<<", "<<Sum_num_trackable<<", "<<Sum_num_golden<<", "<<Sum_num_zero<<std::endl;
  std::cout<<"sumold (2mu): "<<Sum_den_2presel<<", "<<Sum_num_2trackable<<", "<<Sum_num_2golden<<std::endl;
  std::cout<<"sumold (3mu): "<<Sum_den_3presel<<", "<<Sum_num_3trackable<<", "<<Sum_num_3golden<<std::endl;
  std::cout<<"sumold (4mu): "<<Sum_den_4presel<<", "<<Sum_num_4trackable<<", "<<Sum_num_4golden<<std::endl;

  Long64_t den_all,
    den_presel, num_trackable, num_golden,
    num_zero,
    den_2presel, num_2trackable, num_2golden,
    den_3presel, num_3trackable, num_3golden,
    den_4presel, num_4trackable, num_4golden;
  
  read_counts(i,
	      den_all,
	      den_presel, num_trackable, num_golden,
	      num_zero,
	      den_2presel, num_2trackable, num_2golden,
	      den_3presel, num_3trackable, num_3golden,
	      den_4presel, num_4trackable, num_4golden);

  std::cout<<"last (1mu): "<<den_all<<", "<<den_presel<<", "<<num_trackable<<", "<<num_golden<<", "<<num_zero<<std::endl;
  std::cout<<"last (2mu): "<<den_2presel<<", "<<num_2trackable<<", "<<num_2golden<<std::endl;
  std::cout<<"last (3mu): "<<den_3presel<<", "<<num_3trackable<<", "<<num_3golden<<std::endl;
  std::cout<<"last (4mu): "<<den_4presel<<", "<<num_4trackable<<", "<<num_4golden<<std::endl;

  Sum_den_all += den_all;
  Sum_den_presel += den_presel;
  Sum_num_trackable += num_trackable;
  Sum_num_golden += num_golden;
  Sum_num_zero += num_zero;
  Sum_den_2presel += den_2presel;
  Sum_num_2trackable += num_2trackable;
  Sum_num_2golden += num_2golden;
  Sum_den_3presel += den_3presel;
  Sum_num_3trackable += num_3trackable;
  Sum_num_3golden += num_3golden;
  Sum_den_4presel += den_4presel;
  Sum_num_4trackable += num_4trackable;
  Sum_num_4golden += num_4golden;
  //
  std::pair ineff_zero = binomial_eff(Sum_num_zero, Sum_den_all);
  std::pair eff_trackable_all = binomial_eff(Sum_num_trackable, Sum_den_all);
  std::pair eff_trackable_presel = binomial_eff(Sum_num_trackable, Sum_den_presel);
  std::pair eff_golden_presel = binomial_eff(Sum_num_golden, Sum_den_presel);
  std::pair eff_2trackable_presel = binomial_eff(Sum_num_2trackable, Sum_den_2presel);
  std::pair eff_2golden_presel = binomial_eff(Sum_num_2golden, Sum_den_2presel);
  std::pair eff_3trackable_presel = binomial_eff(Sum_num_3trackable, Sum_den_3presel);
  std::pair eff_3golden_presel = binomial_eff(Sum_num_3golden, Sum_den_3presel);
  std::pair eff_4trackable_presel = binomial_eff(Sum_num_4trackable, Sum_den_4presel);
  std::pair eff_4golden_presel = binomial_eff(Sum_num_4golden, Sum_den_4presel);
  //
  
  sum<<"RELATIVE EFFICIENCY Station S_"<<inum
     <<" (in events with Golden Passing Mu in S_"<<iden<<")"<<std::endl;

  sum<< "#all_S"<<iden<<"_passing_golden_events: " << Sum_den_all <<std::endl;
  sum<< "#presel_S"<<iden<<"_passing_golden_events: " << Sum_den_presel <<std::endl;
  sum<< "#S"<<iden<<"_passing_golden_&_S"<<inum<<"_trackable: " << Sum_num_trackable << std::endl;
  sum<< "#passing_golden_in_S"<<iden<<"_&&_S"<<inum<<": " << Sum_num_golden << std::endl;
  sum<< "efficiency S"<<inum<<" trackable / all S"<<iden<<" passing golden = " << eff_trackable_all.first << " +- " << eff_trackable_all.second <<std::endl;
  sum<< "efficiency S"<<inum<<" trackable / presel S"<<iden<<" pass.golden = " << eff_trackable_presel.first << " +- " << eff_trackable_presel.second <<std::endl;
  sum<< "efficiency S"<<inum<<" golden    / presel S"<<iden<<" pass.golden = " << eff_golden_presel.first << " +- " << eff_golden_presel.second <<std::endl;
  sum<< "#S"<<inum<<"_zero: " << Sum_num_zero << std::endl;
  sum<< "fraction ZERO/all : " << ineff_zero.first << " +- " << ineff_zero.second <<std::endl;

  sum<< "#presel_2mu_S"<<iden<<"_passing_golden_events: " << Sum_den_2presel <<std::endl;
  sum<< "#2mu_S"<<iden<<"_passing_golden_&_S"<<inum<<"_trackable: " << Sum_num_2trackable << std::endl;
  sum<< "#2mu_passing_golden_in_S"<<iden<<"_&&_S"<<inum<<": " << Sum_num_2golden << std::endl;
  sum<< "2mu efficiency S"<<inum<<" trackable / presel S"<<iden<<" passing golden = " << eff_2trackable_presel.first << " +- " << eff_2trackable_presel.second <<std::endl;
  sum<< "2mu efficiency S"<<inum<<" golden    / presel S"<<iden<<" pass.golden    = " << eff_2golden_presel.first << " +- " << eff_2golden_presel.second <<std::endl;

  sum<< "#presel_3mu_S"<<iden<<"_passing_golden_events: " << Sum_den_3presel <<std::endl;
  sum<< "#3mu_S"<<iden<<"_passing_golden_&_S"<<inum<<"_trackable: " << Sum_num_3trackable << std::endl;
  sum<< "#3mu_passing_golden_in_S"<<iden<<"_&&_S"<<inum<<": " << Sum_num_3golden << std::endl;
  sum<< "3mu efficiency S"<<inum<<" trackable / presel S"<<iden<<" passing golden = " << eff_3trackable_presel.first << " +- " << eff_3trackable_presel.second <<std::endl;
  sum<< "3mu efficiency S"<<inum<<" golden    / presel S"<<iden<<" pass.golden    = " << eff_3golden_presel.first << " +- " << eff_3golden_presel.second <<std::endl;

  sum<< "#presel_4mu_S"<<iden<<"_passing_golden_events: " << Sum_den_4presel <<std::endl;
  sum<< "#4mu_S"<<iden<<"_passing_golden_&_S"<<inum<<"_trackable: " << Sum_num_4trackable << std::endl;
  sum<< "#4mu_passing_golden_in_S"<<iden<<"_&&_S"<<inum<<": " << Sum_num_4golden << std::endl;
  sum<< "4mu efficiency S"<<inum<<" trackable / presel S"<<iden<<" passing golden = " << eff_4trackable_presel.first << " +- " << eff_4trackable_presel.second <<std::endl;
  sum<< "4mu efficiency S"<<inum<<" golden    / presel S"<<iden<<" pass.golden    = " << eff_4golden_presel.first << " +- " << eff_4golden_presel.second <<std::endl;
    
  std::cout<<"sum (1mu): "<<Sum_den_all<<", "<<Sum_den_presel<<", "<<Sum_num_trackable<<", "<<Sum_num_golden<<", "<<Sum_num_zero<<std::endl;
  std::cout<<"sum (2mu): "<<Sum_den_2presel<<", "<<Sum_num_2trackable<<", "<<Sum_num_2golden<<std::endl;
  std::cout<<"sum (3mu): "<<Sum_den_3presel<<", "<<Sum_num_3trackable<<", "<<Sum_num_3golden<<std::endl;
  std::cout<<"sum (4mu): "<<Sum_den_4presel<<", "<<Sum_num_4trackable<<", "<<Sum_num_4golden<<std::endl;

};

#endif
