#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include "SimStruct.h"
#include <fstream>
#include <cstring>
#include <unordered_set>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

void usage() {
    cerr << "./SimTab-Baseline" << endl
	 << "-d <dataset>" << endl
	 << "-f <filelabel>" << endl
	 << "-algo <algorithm>" << endl
 	 << "[-c <damping factor> (default 0.6)]" << endl
	 << "[-qn <querynum> (default 1)]" << endl
	 << "[-k <k> (default 50)]" << endl;
}

//Check parameters
int check_inc(int i, int max) {
    if (i == max) {
        usage();
        exit(1);
    }
    return i + 1;
}

bool cmp(const pair<int, double>& a, const pair<int, double>& b){ // decreasing order
    return a.second > b.second;
}

int main(int argc, char **argv){
    int i = 1;
    char *endptr;
    string filename;
    string outputname = "-1";
    string filelabel;
    string algo;
    int querynum = -1;
    double eps = 0.1;
    double c = 0.6;
    int k = 50;
    int indexnum=0;
    int isShowResult = 1;
    bool use_index = true;
    bool check_DuplicateEdge = false;
      
    if(argc < 4){
		cout << "argc < 4" << endl;
        usage();
        exit(1);
    }
    while (i < argc) {
        if (!strcmp(argv[i], "-d")) {
            i = check_inc(i, argc);
            filename = argv[i];
        } else if (!strcmp(argv[i], "-f")) {
            i = check_inc(i, argc);
            filelabel = argv[i];
        } else if (!strcmp(argv[i], "-algo")) {
            i = check_inc(i, argc);
            algo = argv[i];
        } else if (!strcmp(argv[i], "-k")) {
			i = check_inc(i, argc);
			k = atoi(argv[i]);
			if (k <= 0 || k >= 5000) {
				cerr << "Invalid k argument" << endl;
                exit(1);
			}
			cout << "k = " << k << endl;
		}
		else if (!strcmp(argv[i], "-c")) {
            i = check_inc(i, argc);
            c = strtod(argv[i], &endptr);
            if ((c == 0 || c > 1) && endptr) {
                cerr << "Invalid c argument" << endl;
                exit(1);
            }
        } else if (!strcmp(argv[i], "-qn")) {
            i = check_inc(i, argc);
            querynum = strtod(argv[i], &endptr);
            if ((querynum < 0) && endptr) {
                cerr << "Invalid querynum argument" << endl;
                exit(1);
            }
        } else if (!strcmp(argv[i], "-o")) {
            i = check_inc(i, argc);
            outputname = argv[i];
        } else {
			cout << "invalid argument: i = " << i << " , argv[i] = " << argv[i] << endl;
            usage();
            exit(1);
        }
        i++;
    }
    if(filename == "" || filelabel == "" || algo == ""){
		cout << "Miss parameters" << endl; 
   		usage();
		exit(1);
    }

    SimStruct sim = SimStruct(filename, filelabel, eps, c);
    
    if(querynum == -1 || querynum > sim.vert)
        querynum = 1;
	
    if(algo == "topk"){
        string queryname = "query/" + filelabel + ".query";
		ifstream fquery;
		fquery.open(queryname.data());
		assert(fquery.is_open());
		cout << "querynum = " << querynum << endl;
		for(int i = 0; i < querynum; i++){
			int nodeId;
            fquery >> nodeId;
            cout << i << ": " << nodeId << endl;
            sim.prefilter(nodeId, k);
			sim.outputCandidate(nodeId, k);			
			stringstream ss_cand_path;
			ss_cand_path << "/home/liuyu/SimTab-release-v1/candidates/" << filelabel << "/k=" << k << "/";
			sim.readCandidatesFromPrefilter(ss_cand_path.str(), nodeId);
			Timer ugape_timer;
			ugape_timer.start();
			sim.UGapE(nodeId, k);
			ugape_timer.end();
			cout << "nodeId = " << nodeId << " , k = " << k << " , ugape_timer.timeCost = " << ugape_timer.timeCost << endl;
        }
		fquery.close();    
		cout << "avg time: " << sim.avg_time / (double) querynum << endl;
		cout << "====done===="<<endl;
    }
    return 0;
}
