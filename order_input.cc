#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <regex>


using namespace std;

////////////////////////////////////////////////////////////////////////////////////
// r1 and r2 are file names coded as:  <start_SuperId>_<end_SuperId>.root
//
bool compare_interval(string r1, string r2) {
  const string numbers("0123456789");

auto pos1 = r1.find('-');
auto pos1_2 = r1.find('.', pos1);  // trova il punto dopo il '-'

unsigned int start_SuperId_1 = std::stoul(r1.substr(pos1 + 1, pos1_2 - pos1 - 1));

auto pos2 = r2.find('-');
auto pos2_2 = r2.find('.', pos1);  // trova il punto dopo il '-'

unsigned int start_SuperId_2 = std::stoul(r2.substr(pos2 + 1, pos2_2 - pos2 - 1));


  if (start_SuperId_1 < start_SuperId_2 ||
      start_SuperId_2 < start_SuperId_1 )
    cout<<"\n ***WARNING: overlapping intervals in chunks "<<r1<<", "<<r2<<endl;
					   
  return (start_SuperId_1 < start_SuperId_2);
}
////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
  
  if(argc != 2) {
    cout<<"ERROR! Usage: "<<argv[0]<<" run_<number>"<<endl;
    return 10;
  }

  string run_nr = argv[1];
  cout<<"\n"<<"INPUT: "<< run_nr <<endl;

  string input = "runs/"+run_nr+"_chunks_dgm.txt";
  ifstream input_file_list(input);
  string i_file;
  int n_files = 0;

  vector<string> files;
  while(getline(input_file_list, i_file)) {
    n_files++;
    //    cout<< i_file << endl;
    files.push_back(i_file);
  }
  cout<<"reading "<<n_files<<" files from: "<<input<<endl;
  cout<<"sorting files according to the SuperID interval" <<endl;
  //=================================================================
  sort(files.begin(), files.end(), compare_interval);
  //=================================================================

  string output = "runs/"+run_nr+"_files.txt";

  ofstream os(output);
  for (int i = 0; i < n_files; i++) {
    os << files.at(i) <<endl;
    //    cout<<i<<". "<<files.at(i)<<endl;
  }
  os.close();
  cout<<" --> "<<output<<endl;
  
  return 0;
}










