#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char *argv[]) {
  
  if(argc != 3) {
    cout<<"ERROR! Usage: "<<argv[0]<<" <run_Nr> <group N files>"<<endl;
    return 10;
  }
  cout<<"\n INPUT: "<< argv[1] <<endl;

  string input_run = argv[1];
  int n_files_per_group = std::stoi(argv[2]);
//  string input = "runs/"+input_run+"_files.txt";
  string input = "runs/"+input_run+"_chunks_dgm.txt";
  ifstream input_file_list(input.c_str());
  string i_file;
  int n_files = 0;
  int n_group = 0;

  string split_name = "inputs/"+input_run+"/inputs_";
  string split_root = split_name + std::to_string(n_group);
  
  ofstream os(split_root.c_str());
  while(getline(input_file_list, i_file)) {
    n_files++;
    cout<< i_file << endl;
    int isplit = n_files%n_files_per_group;
    if (isplit == 0) {
      os.close();
      n_group++;
      split_root = split_name + std::to_string(n_group);
      os.open(split_root.c_str());
    }
    os << i_file <<endl;
  }
  os.close();

  return 0;
}











